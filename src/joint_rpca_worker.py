import argparse
import json
from pathlib import Path

import pandas as pd
from biom import Table
from gemelli.rpca import joint_rpca, feature_correlation_table
from gemelli.utils import filter_ordination  # noqa: F401

from joint_rpca import joint_rpca_scores_to_long


def _feature_table_from_dataframe(df):
    return Table(
        df.to_numpy(dtype=float),
        observation_ids=df.index.astype(str).tolist(),
        sample_ids=df.columns.astype(str).tolist(),
    )


def run_worker(args):
    metabolome_ft = pd.read_csv(args.metabolome, sep="\t", index_col=0)
    omics_ft = pd.read_csv(args.omics, sep="\t", index_col=0)

    table_metabolome = _feature_table_from_dataframe(metabolome_ft.fillna(0))
    table_omics = _feature_table_from_dataframe(omics_ft.fillna(0))

    biplot, distance_matrix, cross_validation_error = joint_rpca(
        [table_metabolome, table_omics],
        min_feature_frequency=args.min_feature_frequency,
        max_iterations=args.max_iterations,
    )

    correlation_matrix = feature_correlation_table(biplot)
    scores, cross_block = joint_rpca_scores_to_long(correlation_matrix)

    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    scores.to_csv(output_dir / "scores.csv", index=False)
    cross_block.to_csv(output_dir / "cross_block.csv")
    correlation_matrix.to_csv(output_dir / "all_feature_scores.csv")
    cross_validation_error.to_csv(output_dir / "cross_validation_error.csv")

    with open(output_dir / "metadata.json", "w") as handle:
        json.dump({"common_samples": metabolome_ft.columns.astype(str).tolist()}, handle)


def main():
    parser = argparse.ArgumentParser(description="Run standalone Gemelli Joint-RPCA for CorrOmics.")
    parser.add_argument("--metabolome", required=True)
    parser.add_argument("--omics", required=True)
    parser.add_argument("--output-dir", required=True)
    parser.add_argument("--max-iterations", type=int, default=5)
    parser.add_argument("--min-feature-frequency", type=float, default=2)
    args = parser.parse_args()
    run_worker(args)


if __name__ == "__main__":
    main()
