# Methods

## Dataset

This project analyzes the immune-cell compartment of the public ulcerative-colitis single-cell atlas from Smillie et al., *Cell* (2019), doi:10.1016/j.cell.2019.06.029. Data are available from the Broad Institute Single Cell Portal under accession SCP259.

The immune-cell data are stored as separate raw files:

- `Imm.barcodes2.tsv`: immune cell barcodes
- `Imm.genes.tsv`: gene names
- `gene_sorted-Imm.matrix.mtx`: sparse gene expression matrix
- `Imm.tsne.txt`: published t-SNE coordinates
- `all.meta2.txt`: cell-level metadata

Raw and processed data files are not tracked in Git because of file size.

## Data Loading

The raw immune-cell expression matrix is loaded from Matrix Market format and combined with immune cell barcodes and gene names. Because the raw matrix is stored as genes by cells, it is transposed before creating a Scanpy `AnnData` object.

Cell-level metadata are matched to immune cell barcodes and stored in `adata.obs`.

## Data Object

The main analysis object is a Scanpy `AnnData` object.

- `adata.X`: gene expression matrix
- `adata.obs`: cell metadata
- `adata.var_names`: gene names
- `adata.obs_names`: cell barcodes

The processed object is saved locally as:

```text
data/processed/immune_cells_with_metadata.h5ad
```

## Donor-aware pseudobulk validation

The original notebooks use cell-level Wilcoxon tests for exploratory marker and disease-associated expression analysis. Because cells from the same biological subject are not independent replicates, confirmatory analysis should operate at the donor level.

The `immune-atlas-pseudobulk` command implements an auditable donor-aware workflow:

1. Select one immune population and two tissue conditions.
2. Validate that the chosen AnnData matrix contains raw, non-negative integer counts.
3. Sum counts within each `Subject` and `Health` combination.
4. Remove donor samples with too few cells and genes with insufficient prevalence.
5. Convert pseudobulk counts to log2 counts per million.
6. Run an exploratory Welch test across donor-level samples and control false discovery rate with Benjamini-Hochberg adjustment.

This test is suitable as a transparent portfolio-scale validation, but a negative-binomial count model such as edgeR or DESeq2 is preferable for a publication-grade analysis. The cell-level and donor-level results must be reported separately.
