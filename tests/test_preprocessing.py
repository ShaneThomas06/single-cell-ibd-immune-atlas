import unittest

import numpy as np
import pandas as pd

from immune_atlas.preprocessing import (
    filter_genes_by_sample_prevalence,
    normalize_total_log1p,
    validate_count_matrix,
)


class PreprocessingTests(unittest.TestCase):
    def test_validate_count_matrix_rejects_log_values(self):
        counts = pd.DataFrame([[0.0, 1.2], [3.0, 4.0]])
        with self.assertRaisesRegex(ValueError, "raw integer counts"):
            validate_count_matrix(counts)

    def test_normalize_total_log1p_preserves_shape_and_equalizes_depth(self):
        counts = pd.DataFrame([[1, 3], [2, 6]], index=["a", "b"], columns=["g1", "g2"])
        result = normalize_total_log1p(counts, target_sum=100)
        reconstructed = np.expm1(result)
        self.assertEqual(result.shape, counts.shape)
        self.assertTrue(np.allclose(reconstructed.sum(axis=1), 100))

    def test_filter_genes_by_sample_prevalence(self):
        counts = pd.DataFrame(
            {"keep": [10, 11, 0], "drop": [10, 0, 0]},
            index=["s1", "s2", "s3"],
        )
        result = filter_genes_by_sample_prevalence(counts, min_count=10, min_samples=2)
        self.assertEqual(result.columns.tolist(), ["keep"])


if __name__ == "__main__":
    unittest.main()
