import importlib.util
import json
import os
import subprocess
import sys
import tempfile
from pathlib import Path

import pandas as pd


JOINT_RPCA_METABOLOME_PREFIX = "metabolome__"
JOINT_RPCA_OMICS_PREFIX = "omics__"


if "CORROMICS_GEMELLI_BACKEND" not in os.environ:
    os.environ["CORROMICS_GEMELLI_BACKEND"] = "disabled" if os.name == "nt" else "local"


def get_gemelli_worker_command():
    """
    Return the local command prefix used to run Python with Gemelli available.
    """

    backend = os.environ.get("CORROMICS_GEMELLI_BACKEND", "disabled").strip().lower()

    if backend == "disabled":
        return None

    if backend != "local":
        raise RuntimeError(
            f"Unsupported CORROMICS_GEMELLI_BACKEND={backend!r}. "
            "Use 'disabled' or 'local'."
        )

    if importlib.util.find_spec("gemelli"):
        return [sys.executable]

    return None

def get_gemelli_worker_python():
    """
    Backward-compatible helper for older code paths.
    Prefer get_gemelli_worker_command() for new code.
    """
    cmd = get_gemelli_worker_command()
    if not cmd:
        return None

    return cmd[0]


def check_joint_rpca_dependencies():
    missing = []

    try:
        worker_cmd = get_gemelli_worker_command()
    except RuntimeError as exc:
        missing.append(str(exc))
        return missing

    if not worker_cmd:
        missing.append(
            "Gemelli worker not found. Install gemelli in this active Linux/WSL environment."
        )

    return missing


def _prefix_index(df, prefix):
    df = df.copy()
    df.index = [f"{prefix}{idx}" for idx in df.index.astype(str)]
    return df


def _strip_prefix(value):
    value = str(value)
    for prefix in [JOINT_RPCA_METABOLOME_PREFIX, JOINT_RPCA_OMICS_PREFIX]:
        if value.startswith(prefix):
            return value[len(prefix):]
    return value


def joint_rpca_scores_to_long(correlation_matrix):
    metabolite_cols = [
        col for col in correlation_matrix.columns.astype(str)
        if col.startswith(JOINT_RPCA_METABOLOME_PREFIX)
    ]
    omics_rows = [
        idx for idx in correlation_matrix.index.astype(str)
        if idx.startswith(JOINT_RPCA_OMICS_PREFIX)
    ]

    cross_block = correlation_matrix.loc[omics_rows, metabolite_cols].copy()
    cross_block.index = [_strip_prefix(idx) for idx in cross_block.index]
    cross_block.columns = [_strip_prefix(col) for col in cross_block.columns]

    long_df = (
        cross_block
        .rename_axis(index="Variable", columns="Feature")
        .stack()
        .reset_index(name="Estimate")
    )

    long_df = long_df[["Feature", "Variable", "Estimate"]]
    long_df["P-value"] = pd.NA
    long_df["BH-Corrected P-Value"] = pd.NA
    long_df["R2"] = pd.NA
    long_df["Method"] = "joint_rpca"

    return long_df, cross_block


def _run_joint_rpca_worker(
    metabolome_for_rpca,
    omics_for_rpca,
    common_samples,
    max_iterations,
    min_feature_frequency,
):
    worker_cmd = get_gemelli_worker_command()
    if not worker_cmd:
        raise ImportError("Gemelli worker Python could not be resolved.")

    worker_script = Path(__file__).with_name("joint_rpca_worker.py")

    with tempfile.TemporaryDirectory(prefix="corromics_joint_rpca_") as tmpdir:
        tmpdir = Path(tmpdir)
        metabolome_path = tmpdir / "metabolome.tsv"
        omics_path = tmpdir / "omics.tsv"
        output_dir = tmpdir / "output"

        metabolome_for_rpca.to_csv(metabolome_path, sep="\t")
        omics_for_rpca.to_csv(omics_path, sep="\t")

        command = worker_cmd + [
            str(worker_script),
            "--metabolome",
            str(metabolome_path),
            "--omics",
            str(omics_path),
            "--output-dir",
            str(output_dir),
            "--max-iterations",
            str(max_iterations),
            "--min-feature-frequency",
            str(min_feature_frequency),
        ]
        completed = subprocess.run(
            command,
            text=True,
            capture_output=True,
            check=False,
        )

        if completed.returncode != 0:
            stderr = completed.stderr.strip()
            stdout = completed.stdout.strip()
            message = stderr or stdout or "Unknown Joint-RPCA worker error."
            raise RuntimeError(message)

        scores = pd.read_csv(output_dir / "scores.csv")
        cross_block = pd.read_csv(output_dir / "cross_block.csv", index_col=0)
        all_feature_scores = pd.read_csv(output_dir / "all_feature_scores.csv", index_col=0)

        metadata_path = output_dir / "metadata.json"
        if metadata_path.exists():
            with open(metadata_path) as handle:
                metadata = json.load(handle)
            common_samples = pd.Index(metadata.get("common_samples", list(common_samples)))

    return {
        "scores": scores,
        "cross_block": cross_block,
        "all_feature_scores": all_feature_scores,
        "common_samples": common_samples,
    }


def run_joint_rpca_for_omics_pair(
    metabolome_ft,
    omics_ft,
    max_iterations=5,
    min_feature_frequency=2,
):
    missing = check_joint_rpca_dependencies()
    if missing:
        raise ImportError(
            "Joint-RPCA requires Gemelli in the active local Linux/WSL environment: "
            + ", ".join(missing)
        )

    metabolome_ft = metabolome_ft.copy()
    omics_ft = omics_ft.copy()
    metabolome_ft.columns = metabolome_ft.columns.astype(str)
    omics_ft.columns = omics_ft.columns.astype(str)

    common_samples = metabolome_ft.columns.intersection(omics_ft.columns)
    if len(common_samples) < 3:
        raise ValueError(f"Not enough shared samples for Joint-RPCA: {len(common_samples)}")

    metabolome_for_rpca = _prefix_index(
        metabolome_ft[common_samples],
        JOINT_RPCA_METABOLOME_PREFIX,
    )
    omics_for_rpca = _prefix_index(
        omics_ft[common_samples],
        JOINT_RPCA_OMICS_PREFIX,
    )

    if (metabolome_for_rpca < 0).any().any() or (omics_for_rpca < 0).any().any():
        raise ValueError("Joint-RPCA input cannot contain negative values.")

    return _run_joint_rpca_worker(
        metabolome_for_rpca=metabolome_for_rpca.fillna(0),
        omics_for_rpca=omics_for_rpca.fillna(0),
        common_samples=common_samples,
        max_iterations=max_iterations,
        min_feature_frequency=min_feature_frequency,
    )
