"""Small, testable preprocessing helpers used by the analysis notebooks."""

from __future__ import annotations

import numpy as np
import pandas as pd


def validate_count_matrix(counts: pd.DataFrame, *, allow_fractional: bool = False) -> None:
    """Validate that a cell/sample-by-gene table can be treated as raw counts."""
    if counts.empty:
        raise ValueError("Count matrix is empty.")
    values = counts.to_numpy(dtype=float, copy=False)
    if not np.isfinite(values).all():
        raise ValueError("Count matrix contains non-finite values.")
    if (values < 0).any():
        raise ValueError("Count matrix contains negative values.")
    if not allow_fractional and not np.allclose(values, np.round(values), atol=1e-6):
        raise ValueError(
            "Count matrix appears normalized or log-transformed; raw integer counts are required."
        )


def normalize_total_log1p(
    counts: pd.DataFrame,
    *,
    target_sum: float = 10_000.0,
) -> pd.DataFrame:
    """Return library-size-normalized log1p expression without changing the input."""
    validate_count_matrix(counts, allow_fractional=True)
    library_sizes = counts.sum(axis=1)
    if (library_sizes <= 0).any():
        raise ValueError("Every cell or sample must have a positive library size.")
    normalized = counts.div(library_sizes, axis=0) * target_sum
    return np.log1p(normalized)


def filter_genes_by_sample_prevalence(
    counts: pd.DataFrame,
    *,
    min_count: float = 10,
    min_samples: int = 2,
) -> pd.DataFrame:
    """Retain genes meeting a minimum count in at least ``min_samples`` samples."""
    if min_count < 0:
        raise ValueError("min_count must be non-negative.")
    if min_samples < 1:
        raise ValueError("min_samples must be at least one.")
    keep = counts.ge(min_count).sum(axis=0) >= min_samples
    return counts.loc[:, keep].copy()
