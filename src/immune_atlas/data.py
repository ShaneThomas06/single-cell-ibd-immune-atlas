from pathlib import Path

import pandas as pd


REQUIRED_RAW_FILES = [
    "Imm.barcodes2.tsv",
    "Imm.genes.tsv",
    "gene_sorted-Imm.matrix.mtx",
    "Imm.tsne.txt",
    "all.meta2.txt",
]


def check_required_raw_files(data_dir: str | Path) -> dict[str, bool]:
    """Return whether each required raw immune-cell data file exists."""
    data_dir = Path(data_dir)
    return {filename: (data_dir / filename).exists() for filename in REQUIRED_RAW_FILES}


def missing_raw_files(data_dir: str | Path) -> list[str]:
    """Return required raw data files that are missing."""
    file_status = check_required_raw_files(data_dir)
    return [filename for filename, exists in file_status.items() if not exists]


def load_processed_adata(path: str | Path):
    """Load a processed AnnData object from an h5ad file."""
    import scanpy as sc

    path = Path(path)

    if not path.exists():
        raise FileNotFoundError(f"Processed AnnData file not found: {path}")

    return sc.read_h5ad(path)


def save_table(df: pd.DataFrame, path: str | Path) -> Path:
    """Save a DataFrame as CSV and create parent directories if needed."""
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(path, index=False)
    return path
