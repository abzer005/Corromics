import importlib.util
import json
import os
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

import pandas as pd


JOINT_RPCA_METABOLOME_PREFIX = "metabolome__"
JOINT_RPCA_OMICS_PREFIX = "omics__"
GEMELLI_WORKER_ENV_NAME = "gemelli-standalone"


def check_joint_rpca_dependencies():
    missing = []
    worker_python = get_gemelli_worker_python()
    if not worker_python:
        missing.append(
            "Gemelli worker Python not found. Install gemelli in this environment, "
            f"create a conda env named {GEMELLI_WORKER_ENV_NAME}, or set GEMELLI_WORKER_PYTHON."
        )
    elif not Path(worker_python).exists():
        missing.append(f"Gemelli worker Python not found: {worker_python}")

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


def _python_from_env_prefix(env_prefix):
    env_prefix = Path(env_prefix).expanduser()
    scripts_dir = "Scripts" if os.name == "nt" else "bin"
    python_name = "python.exe" if os.name == "nt" else "python"
    return env_prefix / scripts_dir / python_name


def _conda_env_prefix_from_cli(env_name):
    for executable_name in ["conda", "mamba", "micromamba"]:
        executable = shutil.which(executable_name)
        if not executable:
            continue

        try:
            completed = subprocess.run(
                [executable, "env", "list", "--json"],
                text=True,
                capture_output=True,
                check=False,
                timeout=15,
            )
        except (OSError, subprocess.TimeoutExpired):
            continue

        if completed.returncode != 0:
            continue

        try:
            envs = json.loads(completed.stdout).get("envs", [])
        except json.JSONDecodeError:
            continue

        for env_prefix in envs:
            env_prefix = Path(env_prefix)
            if env_prefix.name == env_name:
                return env_prefix

    return None


def _candidate_gemelli_worker_pythons(env_name=GEMELLI_WORKER_ENV_NAME):
    env_prefix = _conda_env_prefix_from_cli(env_name)
    if env_prefix:
        yield _python_from_env_prefix(env_prefix)

    for env_root in [
        os.environ.get("CONDA_ENVS_PATH"),
        Path(os.environ["CONDA_PREFIX"]).parent if os.environ.get("CONDA_PREFIX") else None,
        Path.home() / "miniconda3" / "envs",
        Path.home() / "anaconda3" / "envs",
        Path.home() / "mambaforge" / "envs",
        Path.home() / "miniforge3" / "envs",
        Path.home() / "micromamba" / "envs",
        Path("/opt/conda/envs"),
    ]:
        if not env_root:
            continue
        for root in str(env_root).split(os.pathsep):
            yield _python_from_env_prefix(Path(root) / env_name)


def get_gemelli_worker_python():
    configured_python = os.environ.get("GEMELLI_WORKER_PYTHON")
    if configured_python:
        return configured_python

    if importlib.util.find_spec("gemelli"):
        return sys.executable

    for worker_python in _candidate_gemelli_worker_pythons():
        if worker_python.exists():
            return str(worker_python)

    return None


def _run_joint_rpca_worker(
    metabolome_for_rpca,
    omics_for_rpca,
    common_samples,
    max_iterations,
    min_feature_frequency,
):
    worker_python = get_gemelli_worker_python()
    if not worker_python:
        raise ImportError("Gemelli worker Python could not be resolved.")

    worker_script = Path(__file__).with_name("joint_rpca_worker.py")

    with tempfile.TemporaryDirectory(prefix="corromics_joint_rpca_") as tmpdir:
        tmpdir = Path(tmpdir)
        metabolome_path = tmpdir / "metabolome.tsv"
        omics_path = tmpdir / "omics.tsv"
        output_dir = tmpdir / "output"

        metabolome_for_rpca.to_csv(metabolome_path, sep="\t")
        omics_for_rpca.to_csv(omics_path, sep="\t")

        command = [
            worker_python,
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
            "Joint-RPCA requires the optional standalone Gemelli worker environment: "
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
