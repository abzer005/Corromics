import time
import numpy as np
import pandas as pd


def distance_correlation_1d(x, y):
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)

    if x.shape[0] != y.shape[0]:
        raise ValueError("x and y must have the same length.")

    n = x.shape[0]
    if n < 3:
        return np.nan

    a = np.abs(x[:, None] - x[None, :])
    b = np.abs(y[:, None] - y[None, :])

    A = a - a.mean(axis=0)[None, :] - a.mean(axis=1)[:, None] + a.mean()
    B = b - b.mean(axis=0)[None, :] - b.mean(axis=1)[:, None] + b.mean()

    dcov2 = np.mean(A * B)
    dvarx = np.mean(A * A)
    dvary = np.mean(B * B)

    if dvarx <= 0 or dvary <= 0:
        return 0.0

    return np.sqrt(dcov2) / np.sqrt(np.sqrt(dvarx * dvary))


def compute_distance_correlation(
    metabolome_ft: pd.DataFrame,
    omics_ft: pd.DataFrame,
    progress_callback=None,
):
    common_samples = sorted(set(metabolome_ft.columns) & set(omics_ft.columns))

    if len(common_samples) < 3:
        raise ValueError("At least 3 shared samples are required for distance correlation.")

    metabolome_sub = metabolome_ft[common_samples].copy()
    omics_sub = omics_ft[common_samples].copy()

    total_pairs = metabolome_sub.shape[0] * omics_sub.shape[0]
    completed_pairs = 0
    start_time = time.time()

    results = []

    for met_id, met_row in metabolome_sub.iterrows():
        x = met_row.to_numpy(dtype=float)

        if np.all(~np.isfinite(x)) or np.nanstd(x) == 0:
            completed_pairs += omics_sub.shape[0]
            if progress_callback is not None:
                progress_callback(completed_pairs, total_pairs, start_time)
            continue

        for omics_id, omics_row in omics_sub.iterrows():
            y = omics_row.to_numpy(dtype=float)

            if not (np.all(~np.isfinite(y)) or np.nanstd(y) == 0):
                mask = np.isfinite(x) & np.isfinite(y)
                if mask.sum() >= 3:
                    score = distance_correlation_1d(x[mask], y[mask])

                    results.append({
                        "Feature": met_id,
                        "Variable": omics_id,
                        "Estimate": float(score),
                        "P-value": np.nan,
                        "BH-Corrected P-Value": np.nan,
                        "Method": "distance_correlation",
                        "R2": np.nan,
                    })

            completed_pairs += 1

            if progress_callback is not None and (
                completed_pairs % 100 == 0 or completed_pairs == total_pairs
            ):
                progress_callback(completed_pairs, total_pairs, start_time)

    results_df = pd.DataFrame(results)
    return results_df, common_samples


def run_distance_correlation(
    metabolome_ft: pd.DataFrame,
    omics_ft: pd.DataFrame,
    progress_callback=None,
):
    results_df, common_samples = compute_distance_correlation(
        metabolome_ft=metabolome_ft,
        omics_ft=omics_ft,
        progress_callback=progress_callback,
    )
    return results_df, common_samples