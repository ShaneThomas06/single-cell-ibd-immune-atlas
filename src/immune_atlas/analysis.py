"""Reusable composition and donor-aware pseudobulk analysis helpers."""

from __future__ import annotations

import numpy as np
import pandas as pd
from scipy.stats import ttest_ind

from immune_atlas.preprocessing import validate_count_matrix


def composition_by_condition(
    metadata: pd.DataFrame,
    *,
    cluster_col: str = "Cluster",
    condition_col: str = "Health",
) -> pd.DataFrame:
    """Summarize cell counts and within-condition percentages by cluster."""
    missing = {cluster_col, condition_col}.difference(metadata.columns)
    if missing:
        raise ValueError(f"Metadata is missing columns: {', '.join(sorted(missing))}")
    counts = (
        metadata.groupby([condition_col, cluster_col], observed=True)
        .size()
        .rename("n_cells")
        .reset_index()
    )
    totals = counts.groupby(condition_col)["n_cells"].transform("sum")
    counts["percentage"] = counts["n_cells"] / totals * 100
    return counts.sort_values([condition_col, "n_cells"], ascending=[True, False])


def aggregate_pseudobulk_counts(
    counts: pd.DataFrame,
    metadata: pd.DataFrame,
    *,
    donor_col: str = "Subject",
    condition_col: str = "Health",
    min_cells: int = 10,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Sum raw single-cell counts into donor-by-condition pseudobulk samples."""
    if min_cells < 1:
        raise ValueError("min_cells must be at least one.")
    missing_obs = counts.index.difference(metadata.index)
    if not missing_obs.empty:
        raise ValueError("Metadata is missing one or more count-matrix rows.")
    metadata = metadata.loc[counts.index]
    missing_cols = {donor_col, condition_col}.difference(metadata.columns)
    if missing_cols:
        raise ValueError(f"Metadata is missing columns: {', '.join(sorted(missing_cols))}")
    validate_count_matrix(counts)

    groupers = [metadata[donor_col], metadata[condition_col]]
    aggregated = counts.groupby(groupers, observed=True).sum()
    cell_counts = metadata.groupby(groupers, observed=True).size().rename("n_cells")
    keep = cell_counts >= min_cells
    aggregated = aggregated.loc[keep[keep].index]
    cell_counts = cell_counts.loc[aggregated.index]
    if aggregated.empty:
        raise ValueError("No pseudobulk samples remain after the minimum-cell filter.")

    sample_metadata = aggregated.index.to_frame(index=False)
    sample_metadata["sample_id"] = (
        sample_metadata[donor_col].astype(str)
        + "__"
        + sample_metadata[condition_col].astype(str)
    )
    sample_metadata["n_cells"] = cell_counts.to_numpy()
    sample_metadata = sample_metadata.set_index("sample_id")
    aggregated.index = sample_metadata.index
    return aggregated, sample_metadata


def log2_cpm(counts: pd.DataFrame, *, prior_count: float = 0.5) -> pd.DataFrame:
    """Convert pseudobulk counts to log2 counts per million."""
    validate_count_matrix(counts)
    if prior_count <= 0:
        raise ValueError("prior_count must be positive.")
    library_sizes = counts.sum(axis=1)
    if (library_sizes <= 0).any():
        raise ValueError("Every pseudobulk sample must have a positive library size.")
    cpm = counts.div(library_sizes, axis=0) * 1_000_000
    return np.log2(cpm + prior_count)


def benjamini_hochberg(pvalues: np.ndarray | pd.Series) -> np.ndarray:
    """Adjust p-values with the Benjamini-Hochberg false-discovery-rate method."""
    values = np.asarray(pvalues, dtype=float)
    adjusted = np.full(values.shape, np.nan, dtype=float)
    valid = np.isfinite(values)
    if not valid.any():
        return adjusted
    valid_values = values[valid]
    order = np.argsort(valid_values)
    ranked = valid_values[order]
    n_tests = ranked.size
    ranked_adjusted = ranked * n_tests / np.arange(1, n_tests + 1)
    ranked_adjusted = np.minimum.accumulate(ranked_adjusted[::-1])[::-1]
    restored = np.empty_like(ranked_adjusted)
    restored[order] = np.clip(ranked_adjusted, 0, 1)
    adjusted[valid] = restored
    return adjusted


def welch_pseudobulk_de(
    expression: pd.DataFrame,
    sample_metadata: pd.DataFrame,
    *,
    test_group: str,
    reference_group: str,
    condition_col: str = "Health",
    min_replicates: int = 2,
) -> pd.DataFrame:
    """Run an exploratory donor-level Welch test on log2-CPM pseudobulk values."""
    if condition_col not in sample_metadata.columns:
        raise ValueError(f"Sample metadata is missing column: {condition_col}")
    missing_samples = expression.index.difference(sample_metadata.index)
    if not missing_samples.empty:
        raise ValueError("Sample metadata is missing one or more expression rows.")
    sample_metadata = sample_metadata.loc[expression.index]
    test_mask = sample_metadata[condition_col].astype(str).eq(test_group)
    reference_mask = sample_metadata[condition_col].astype(str).eq(reference_group)
    if test_mask.sum() < min_replicates or reference_mask.sum() < min_replicates:
        raise ValueError(
            f"Need at least {min_replicates} donor-level samples in both comparison groups."
        )

    test_values = expression.loc[test_mask].to_numpy(dtype=float)
    reference_values = expression.loc[reference_mask].to_numpy(dtype=float)
    statistic, pvalue = ttest_ind(
        test_values,
        reference_values,
        axis=0,
        equal_var=False,
        nan_policy="omit",
    )
    test_mean = np.nanmean(test_values, axis=0)
    reference_mean = np.nanmean(reference_values, axis=0)
    pvalue = np.where(np.isfinite(pvalue), pvalue, 1.0)
    result = pd.DataFrame(
        {
            "gene": expression.columns,
            f"mean_{test_group}": test_mean,
            f"mean_{reference_group}": reference_mean,
            "log2_fold_change": test_mean - reference_mean,
            "statistic": statistic,
            "pvalue": pvalue,
            "pvalue_adj": benjamini_hochberg(pvalue),
        }
    )
    return result.sort_values(["pvalue_adj", "pvalue", "gene"]).reset_index(drop=True)
