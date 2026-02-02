"""Reference data loading utilities."""

import os
import sys
from pathlib import Path

import pandas as pd


def get_data_path() -> Path:
    """
    Get the path to OptiType reference data.

    The data can be located in:
    1. Package installation directory (share/optitype/data)
    2. Development directory (data/)
    3. Custom path via OPTITYPE_DATA environment variable

    Returns:
        Path to the data directory.

    Raises:
        FileNotFoundError: If no valid data directory is found.
    """
    # Check environment variable first
    if env_path := os.environ.get("OPTITYPE_DATA"):
        data_path = Path(env_path)
        if data_path.exists():
            return data_path

    # Check for installed package data
    # When installed via pip, data goes to share/optitype/data
    if sys.prefix != sys.base_prefix:  # In a virtual environment
        share_path = Path(sys.prefix) / "share" / "optitype" / "data"
        if share_path.exists():
            return share_path

    # Also check site-packages share location
    for site_packages in sys.path:
        share_path = Path(site_packages).parent / "share" / "optitype" / "data"
        if share_path.exists():
            return share_path

    # Development mode: check relative to package
    package_dir = Path(__file__).parent.parent.parent.parent
    dev_data_path = package_dir / "data"
    if dev_data_path.exists():
        return dev_data_path

    # Check for data next to the package source
    src_data_path = Path(__file__).parent.parent.parent.parent / "data"
    if src_data_path.exists():
        return src_data_path

    raise FileNotFoundError(
        "Could not find OptiType reference data. "
        "Set OPTITYPE_DATA environment variable to the data directory path, "
        "or ensure the package is properly installed."
    )


def load_reference_data() -> tuple[pd.DataFrame, pd.DataFrame]:
    """
    Load the allele reference data from Parquet files.

    Returns:
        Tuple of (table, features) DataFrames containing allele information.

    Raises:
        FileNotFoundError: If the data files are not found.
    """
    data_path = get_data_path()
    alleles_path = data_path / "alleles"

    if not alleles_path.exists():
        raise FileNotFoundError(
            f"Allele data directory not found at {alleles_path}. "
            "Ensure the data has been converted to Parquet format."
        )

    table_path = alleles_path / "table.parquet"
    features_path = alleles_path / "features.parquet"

    if not table_path.exists() or not features_path.exists():
        raise FileNotFoundError(
            f"Parquet files not found in {alleles_path}. "
            "Expected table.parquet and features.parquet."
        )

    table = pd.read_parquet(table_path)
    features = pd.read_parquet(features_path)

    return table, features


def get_reference_fasta(seq_type: str) -> Path:
    """
    Get the path to the HLA reference FASTA file.

    Args:
        seq_type: Either 'dna' or 'rna'.

    Returns:
        Path to the reference FASTA file.

    Raises:
        ValueError: If seq_type is not 'dna' or 'rna'.
        FileNotFoundError: If the reference file is not found.
    """
    if seq_type not in ("dna", "rna"):
        raise ValueError(f"seq_type must be 'dna' or 'rna', got {seq_type!r}")

    data_path = get_data_path()
    fasta_name = f"hla_reference_{seq_type}.fasta"
    fasta_path = data_path / fasta_name

    if not fasta_path.exists():
        raise FileNotFoundError(f"Reference FASTA not found: {fasta_path}")

    return fasta_path
