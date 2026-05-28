import importlib.util
import sys
import time
from concurrent.futures import ThreadPoolExecutor
from functools import lru_cache
from pathlib import Path

import numpy as np
import pandas as pd
import psutil

SPARCC_RESTRICTED_MAX_FEATURES = 1500
SPARCC_LOCAL_WARNING_FEATURES = 2500
SPARCC_BASE_OVERHEAD_FACTOR = 20
SPARCC_ITERATION_OVERHEAD_FACTOR = 5

def get_sparcc_script_path() -> Path:
    """
    Return absolute path to the vendored SparCC3 runtime packaged with the app.
    Works locally + when deployed.
    """
    here = Path(__file__).resolve().parent
    sparcc_path = here / "vendor" / "sparcc3" / "SparCC.py"

    if not sparcc_path.exists():
        raise FileNotFoundError(f"SparCC script not found at: {sparcc_path}")

    return sparcc_path


def assess_sparcc_run_size(metabolome_ft, omics_ft, is_restricted_mode, iter_num=1):
    """
    Estimate SparCC matrix size and decide whether to warn or block the run.
    """
    total_features = metabolome_ft.shape[0] + omics_ft.shape[0]
    matrix_gb = (total_features ** 2 * 8) / (1024 ** 3)
    available_memory_gb = psutil.virtual_memory().available / (1024 ** 3)
    memory_pressure_gb = matrix_gb * (
        SPARCC_BASE_OVERHEAD_FACTOR + (iter_num * SPARCC_ITERATION_OVERHEAD_FACTOR)
    )
    memory_pressure_ratio = (
        memory_pressure_gb / available_memory_gb
        if available_memory_gb > 0
        else float("inf")
    )

    if (
        is_restricted_mode and total_features > SPARCC_RESTRICTED_MAX_FEATURES
        or memory_pressure_ratio >= 0.70
        or (total_features > SPARCC_LOCAL_WARNING_FEATURES and iter_num >= 3)
    ):
        risk = "High"
    elif (
        memory_pressure_ratio >= 0.35
        or total_features > SPARCC_LOCAL_WARNING_FEATURES
        or iter_num >= 4
    ):
        risk = "Medium"
    else:
        risk = "Low"

    return {
        "total_features": total_features,
        "matrix_gb": matrix_gb,
        "available_memory_gb": available_memory_gb,
        "memory_pressure_gb": memory_pressure_gb,
        "memory_pressure_ratio": memory_pressure_ratio,
        "risk": risk,
        "blocked": is_restricted_mode and total_features > SPARCC_RESTRICTED_MAX_FEATURES,
        "warn": total_features > SPARCC_LOCAL_WARNING_FEATURES,
        "restricted_limit": SPARCC_RESTRICTED_MAX_FEATURES,
        "local_warning_limit": SPARCC_LOCAL_WARNING_FEATURES,
    }

@lru_cache(maxsize=1)
def load_sparcc3_module(sparcc3_script_path: str | Path):
    sparcc3_script_path = Path(sparcc3_script_path).resolve()
    sparcc3_dir = str(sparcc3_script_path.parent)
    if sparcc3_dir not in sys.path:
        sys.path.insert(0, sparcc3_dir)

    spec = importlib.util.spec_from_file_location(
        "corromics_bundled_sparcc3",
        sparcc3_script_path,
    )
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module

def compute_sparcc3_corr(
    X_features_by_samples: pd.DataFrame,
    sparcc3_script_path: str | Path,
    iter_num: int = 1,
    progress_callback=None,
    estimated_seconds: int | float = 60,
) -> pd.DataFrame:
    """
    Run SparCC3 and return the feature×feature correlation matrix.

    Parameters
    ----------
    X_features_by_samples : pd.DataFrame
        Input table with rows = features and columns = samples.
    sparcc3_script_path : str | Path
        Path to SparCC.py

    Returns
    -------
    pd.DataFrame
        SparCC correlation matrix.
    """
    sparcc3_script_path = Path(sparcc3_script_path)

    if not sparcc3_script_path.exists():
        raise FileNotFoundError(f"SparCC script not found: {sparcc3_script_path}")

    if X_features_by_samples.isnull().values.any():
        raise ValueError("Input contains NaNs.")

    if (X_features_by_samples < 0).values.any():
        raise ValueError("Input contains negative values.")

    sparcc_input = X_features_by_samples.copy()
    sparcc3_module = load_sparcc3_module(sparcc3_script_path)
    counts_samples_by_features = sparcc_input.T
    feature_labels = sparcc_input.index.astype(str)

    def _run_sparcc():
        corr_values, _ = sparcc3_module.main(
            counts_samples_by_features,
            method="SparCC",
            iter=iter_num,
            oprint=False,
        )
        return pd.DataFrame(corr_values, index=feature_labels, columns=feature_labels)

    start_time = time.time()
    with ThreadPoolExecutor(max_workers=1) as executor:
        future = executor.submit(_run_sparcc)
        while not future.done():
            if progress_callback is not None:
                progress_callback(start_time, estimated_seconds)
            time.sleep(0.5)

        corr = future.result()

    return corr

def sparcc_corr_to_results_dict(
    corr: pd.DataFrame,
    metabolome_names: pd.Index,
    asv_names: pd.Index,
) -> dict[str, pd.DataFrame]:
    """
    Convert SparCC feature×feature corr matrix into CorrOmics output format:
    {asv: DataFrame(index=metabolites, columns=[Estimate, P-value, BH..., R2])}
    """
    # Ensure all requested names exist in corr
    missing_mets = [m for m in metabolome_names if m not in corr.index]
    missing_asvs = [a for a in asv_names if a not in corr.index]
    if missing_mets or missing_asvs:
        raise ValueError(
            f"SparCC corr matrix missing features. "
            f"Missing metabolites: {missing_mets[:5]}{'...' if len(missing_mets)>5 else ''}; "
            f"Missing ASVs: {missing_asvs[:5]}{'...' if len(missing_asvs)>5 else ''}"
        )

    results = {}
    for asv in asv_names:
        est = corr.loc[metabolome_names, asv].astype(float).values
        out = pd.DataFrame(
            {
                "Estimate": est,
                "P-value": np.nan,
                "BH-Corrected P-Value": np.nan,
                "R2": np.nan,
            },
            index=metabolome_names,
        )
        results[str(asv)] = out
    return results

def run_sparcc_for_omics_pair(
    metabolome_ft: pd.DataFrame,
    omics_ft: pd.DataFrame,
    sparcc_script_path: str | Path,
    iter_num: int = 1,
    progress_callback=None,
    estimated_seconds: int | float = 60,
):
    metabolome_ft = metabolome_ft.copy()
    omics_ft = omics_ft.copy()

    metabolome_ft.columns = metabolome_ft.columns.astype(str)
    omics_ft.columns = omics_ft.columns.astype(str)
    metabolome_ft.index = metabolome_ft.index.astype(str)
    omics_ft.index = omics_ft.index.astype(str)

    common_samples = metabolome_ft.columns.intersection(omics_ft.columns)
    if len(common_samples) < 3:
        raise ValueError(f"Not enough shared samples for SparCC: {len(common_samples)}")

    combined_for_sparcc = pd.concat(
        [
            metabolome_ft[common_samples],
            omics_ft[common_samples],
        ],
        axis=0,
    )

    corr = compute_sparcc3_corr(
        X_features_by_samples=combined_for_sparcc,
        sparcc3_script_path=sparcc_script_path,
        iter_num=iter_num,
        progress_callback=progress_callback,
        estimated_seconds=estimated_seconds,
    )

    corr.index = corr.index.astype(str)
    corr.columns = corr.columns.astype(str)

    results = sparcc_corr_to_results_dict(
        corr=corr,
        metabolome_names=metabolome_ft.index.astype(str),
        asv_names=omics_ft.index.astype(str),
    )

    return corr, results, common_samples
