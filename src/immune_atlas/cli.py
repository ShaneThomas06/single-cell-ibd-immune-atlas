"""Command-line entry points for reproducible analyses."""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import sparse

from immune_atlas.analysis import log2_cpm, welch_pseudobulk_de
from immune_atlas.preprocessing import filter_genes_by_sample_prevalence


def _aggregate_anndata(
    adata,
    *,
    cell_type: str,
    test_group: str,
    reference_group: str,
    cell_type_col: str,
    donor_col: str,
    condition_col: str,
    counts_layer: str,
    min_cells: int,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    required = {cell_type_col, donor_col, condition_col}
    missing = required.difference(adata.obs.columns)
    if missing:
        raise ValueError(f"AnnData metadata is missing columns: {', '.join(sorted(missing))}")
    mask = adata.obs[cell_type_col].astype(str).eq(cell_type) & adata.obs[
        condition_col
    ].astype(str).isin([test_group, reference_group])
    subset = adata[mask]
    if subset.n_obs == 0:
        raise ValueError("No cells match the requested cell type and conditions.")

    if counts_layer == "X":
        matrix = subset.X
    elif counts_layer in subset.layers:
        matrix = subset.layers[counts_layer]
    else:
        raise ValueError(f"Counts layer '{counts_layer}' was not found. Use --counts-layer X if X is raw.")
    observed = matrix.data if sparse.issparse(matrix) else np.asarray(matrix)
    if observed.size and (
        not np.isfinite(observed).all()
        or (observed < 0).any()
        or not np.allclose(observed, np.round(observed), atol=1e-6)
    ):
        raise ValueError("Selected matrix is not raw, non-negative integer count data.")

    group_frame = subset.obs[[donor_col, condition_col]].astype(str)
    group_codes, unique_groups = pd.factorize(pd.MultiIndex.from_frame(group_frame), sort=True)
    group_index = pd.MultiIndex.from_tuples(
        unique_groups.tolist(),
        names=[donor_col, condition_col],
    )
    indicator = sparse.csr_matrix(
        (np.ones(subset.n_obs), (group_codes, np.arange(subset.n_obs))),
        shape=(len(group_index), subset.n_obs),
    )
    matrix = sparse.csr_matrix(matrix)
    aggregated = indicator @ matrix
    cell_counts = np.bincount(group_codes, minlength=len(group_index))
    keep = cell_counts >= min_cells
    if not keep.any():
        raise ValueError("No donor-level samples remain after the minimum-cell filter.")

    group_index = group_index[keep]
    sample_metadata = group_index.to_frame(index=False)
    sample_metadata["sample_id"] = (
        sample_metadata[donor_col] + "__" + sample_metadata[condition_col]
    )
    sample_metadata["n_cells"] = cell_counts[keep]
    sample_metadata = sample_metadata.set_index("sample_id")
    counts = pd.DataFrame(
        aggregated[keep].toarray(),
        index=sample_metadata.index,
        columns=subset.var_names.astype(str),
    )
    return counts, sample_metadata


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Run exploratory donor-aware pseudobulk differential expression."
    )
    parser.add_argument("--input", required=True, help="Processed AnnData .h5ad file")
    parser.add_argument("--output-dir", required=True)
    parser.add_argument("--cell-type", default="Macrophages")
    parser.add_argument("--test-group", default="Inflamed")
    parser.add_argument("--reference-group", default="Healthy")
    parser.add_argument("--cell-type-col", default="Cluster")
    parser.add_argument("--donor-col", default="Subject")
    parser.add_argument("--condition-col", default="Health")
    parser.add_argument("--counts-layer", default="X")
    parser.add_argument("--min-cells", type=int, default=20)
    parser.add_argument("--min-count", type=float, default=10)
    parser.add_argument("--min-samples", type=int, default=2)
    return parser


def main() -> None:
    args = build_parser().parse_args()
    import scanpy as sc

    adata = sc.read_h5ad(args.input)
    counts, sample_metadata = _aggregate_anndata(
        adata,
        cell_type=args.cell_type,
        test_group=args.test_group,
        reference_group=args.reference_group,
        cell_type_col=args.cell_type_col,
        donor_col=args.donor_col,
        condition_col=args.condition_col,
        counts_layer=args.counts_layer,
        min_cells=args.min_cells,
    )
    counts = filter_genes_by_sample_prevalence(
        counts,
        min_count=args.min_count,
        min_samples=args.min_samples,
    )
    expression = log2_cpm(counts)
    results = welch_pseudobulk_de(
        expression,
        sample_metadata,
        test_group=args.test_group,
        reference_group=args.reference_group,
        condition_col=args.condition_col,
    )

    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    counts.to_csv(output_dir / "pseudobulk_counts.csv")
    sample_metadata.to_csv(output_dir / "pseudobulk_samples.csv")
    results.to_csv(output_dir / "pseudobulk_differential_expression.csv", index=False)
    print(f"Saved donor-aware pseudobulk outputs to {output_dir}")


if __name__ == "__main__":
    main()
