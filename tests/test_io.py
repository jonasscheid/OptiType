"""Tests for I/O utilities."""

import pytest
import pandas as pd

from optitype.io.data import get_data_path, load_reference_data, get_reference_fasta


def test_get_data_path():
    """Test that data path can be found."""
    data_path = get_data_path()
    assert data_path.exists()
    assert data_path.is_dir()


def test_load_reference_data():
    """Test loading reference data from Parquet files."""
    table, features = load_reference_data()

    assert isinstance(table, pd.DataFrame)
    assert isinstance(features, pd.DataFrame)

    # Check table structure
    assert "id" in table.columns
    assert "4digit" in table.columns
    assert "locus" in table.columns
    assert len(table) > 0

    # Check features structure
    assert "id" in features.columns
    assert "feature" in features.columns
    assert "number" in features.columns
    assert len(features) > 0


def test_get_reference_fasta():
    """Test getting reference FASTA paths."""
    dna_ref = get_reference_fasta("dna")
    assert dna_ref.exists()
    assert "dna" in dna_ref.name

    rna_ref = get_reference_fasta("rna")
    assert rna_ref.exists()
    assert "rna" in rna_ref.name


def test_get_reference_fasta_invalid():
    """Test that invalid seq_type raises ValueError."""
    with pytest.raises(ValueError):
        get_reference_fasta("invalid")
