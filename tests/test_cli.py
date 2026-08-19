import unittest

import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse

from immune_atlas.cli import _aggregate_anndata


class CliAggregationTests(unittest.TestCase):
    def test_sparse_anndata_is_aggregated_without_cell_level_densification(self):
        matrix = sparse.csr_matrix(
            np.array(
                [
                    [1, 8],
                    [2, 7],
                    [3, 6],
                    [4, 5],
                    [5, 4],
                    [6, 3],
                    [7, 2],
                    [8, 1],
                ]
            )
        )
        obs = pd.DataFrame(
            {
                "Cluster": ["Macrophages"] * 8,
                "Subject": ["H1", "H1", "H2", "H2", "I1", "I1", "I2", "I2"],
                "Health": ["Healthy"] * 4 + ["Inflamed"] * 4,
            },
            index=[f"c{i}" for i in range(8)],
        )
        adata = ad.AnnData(X=matrix, obs=obs, var=pd.DataFrame(index=["g1", "g2"]))
        counts, sample_metadata = _aggregate_anndata(
            adata,
            cell_type="Macrophages",
            test_group="Inflamed",
            reference_group="Healthy",
            cell_type_col="Cluster",
            donor_col="Subject",
            condition_col="Health",
            counts_layer="X",
            min_cells=2,
        )
        self.assertEqual(counts.loc["H1__Healthy", "g1"], 3)
        self.assertEqual(counts.loc["I2__Inflamed", "g2"], 3)
        self.assertTrue(sample_metadata["n_cells"].eq(2).all())


if __name__ == "__main__":
    unittest.main()
