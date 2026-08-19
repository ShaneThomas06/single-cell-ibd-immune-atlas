import unittest

import numpy as np
import pandas as pd

from immune_atlas.analysis import (
    aggregate_pseudobulk_counts,
    benjamini_hochberg,
    composition_by_condition,
    log2_cpm,
    welch_pseudobulk_de,
)


class AnalysisTests(unittest.TestCase):
    def test_composition_by_condition_returns_within_condition_percentages(self):
        metadata = pd.DataFrame(
            {
                "Health": ["Healthy", "Healthy", "Inflamed"],
                "Cluster": ["T", "B", "T"],
            }
        )
        result = composition_by_condition(metadata)
        totals = result.groupby("Health")["percentage"].sum()
        self.assertTrue(np.allclose(totals, 100))

    def test_aggregate_pseudobulk_counts_sums_cells_by_donor_and_condition(self):
        cells = [f"c{i}" for i in range(8)]
        counts = pd.DataFrame(
            {"g1": [1, 2, 3, 4, 5, 6, 7, 8], "g2": [8, 7, 6, 5, 4, 3, 2, 1]},
            index=cells,
        )
        metadata = pd.DataFrame(
            {
                "Subject": ["H1", "H1", "H2", "H2", "I1", "I1", "I2", "I2"],
                "Health": ["Healthy"] * 4 + ["Inflamed"] * 4,
            },
            index=cells,
        )
        pseudobulk, sample_metadata = aggregate_pseudobulk_counts(
            counts,
            metadata,
            min_cells=2,
        )
        self.assertEqual(pseudobulk.loc["H1__Healthy", "g1"], 3)
        self.assertEqual(pseudobulk.loc["I2__Inflamed", "g2"], 3)
        self.assertTrue(sample_metadata["n_cells"].eq(2).all())

    def test_benjamini_hochberg_known_values(self):
        adjusted = benjamini_hochberg(np.array([0.01, 0.04, 0.03, 0.002]))
        self.assertTrue(np.allclose(adjusted, [0.02, 0.04, 0.04, 0.008]))

    def test_welch_pseudobulk_de_prioritizes_shifted_gene(self):
        expression = pd.DataFrame(
            {
                "shifted": [2.0, 2.1, 1.9, 8.0, 8.2, 7.8],
                "stable": [5.0, 5.1, 4.9, 5.0, 5.1, 4.9],
            },
            index=["h1", "h2", "h3", "i1", "i2", "i3"],
        )
        metadata = pd.DataFrame(
            {"Health": ["Healthy"] * 3 + ["Inflamed"] * 3},
            index=expression.index,
        )
        result = welch_pseudobulk_de(
            expression,
            metadata,
            test_group="Inflamed",
            reference_group="Healthy",
        )
        shifted = result.set_index("gene").loc["shifted"]
        self.assertGreater(shifted["log2_fold_change"], 5)
        self.assertEqual(result.iloc[0]["gene"], "shifted")

    def test_log2_cpm_returns_finite_values(self):
        counts = pd.DataFrame({"g1": [10, 20], "g2": [30, 40]}, index=["s1", "s2"])
        result = log2_cpm(counts)
        self.assertTrue(np.isfinite(result.to_numpy()).all())


if __name__ == "__main__":
    unittest.main()
