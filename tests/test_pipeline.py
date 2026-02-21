"""Smoke tests for the pipeline and API."""

import shutil
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from optitype import hlatyper as ht
from optitype.io.data import load_reference_data
from optitype.pipeline import PipelineConfig, run_pipeline, load_config, _is_frequent, _get_4digit, _get_types

# Check for external dependencies
HAS_RAZERS3 = shutil.which("razers3") is not None
HAS_SOLVER = shutil.which("glpsol") is not None or shutil.which("cbc") is not None

# Test data paths
TEST_RNA_1 = Path(__file__).parent.parent / "test" / "rna" / "CRC_81_N_1_fished.fastq"
TEST_RNA_2 = Path(__file__).parent.parent / "test" / "rna" / "CRC_81_N_2_fished.fastq"


class TestPipelineConfig:
    """Tests for PipelineConfig and load_config."""

    def test_default_config(self):
        config = PipelineConfig()
        assert config.razers3_path == "razers3"
        assert config.solver == "glpk"
        assert config.mapping_threads == 4
        assert config.ilp_threads == 1
        assert config.delete_bam is True

    def test_load_config_no_file(self):
        config = load_config(None)
        assert config.razers3_path == "razers3"
        assert config.solver == "glpk"

    def test_load_config_missing_file(self):
        config = load_config("/nonexistent/config.ini")
        assert config.razers3_path == "razers3"


class TestHelperFunctions:
    """Tests for pipeline helper functions."""

    def test_is_frequent(self):
        table, _ = load_reference_data()
        # Find a known frequent allele in the table
        for idx in table.index:
            if table.loc[idx]["4digit"] == "A*02:01" and table.loc[idx]["flags"] == 0:
                assert _is_frequent(idx, table)
                break

    def test_get_4digit(self):
        table, _ = load_reference_data()
        # Pick any allele and check it returns a string with '*'
        idx = table.index[0]
        result = _get_4digit(idx, table)
        assert isinstance(result, str)
        assert "*" in result

    def test_get_types_with_string(self):
        table, _ = load_reference_data()
        idx = table.index[0]
        result = _get_types(idx, table)
        assert isinstance(result, str)

    def test_get_types_with_non_string(self):
        table, _ = load_reference_data()
        assert _get_types(3.14, table) == 3.14
        assert np.isnan(_get_types(float("nan"), table))


class TestHlatyper:
    """Tests for hlatyper matrix operations."""

    def test_get_compact_model(self):
        df = pd.DataFrame(
            [[1, 0, 1], [1, 0, 1], [0, 1, 0], [0, 0, 0]],
            columns=["A", "B", "C"],
            index=["r1", "r2", "r3", "r4"],
        )
        unique_mtx, occ = ht.get_compact_model(df)
        # r4 is all-zero, removed. r1 and r2 are duplicates.
        assert unique_mtx.shape[0] == 2
        assert sum(occ.values()) == 3  # 2 for r1/r2 + 1 for r3

    def test_mtx_to_sparse_dict(self):
        df = pd.DataFrame(
            [[1, 0], [0, 1]],
            columns=["A", "B"],
            index=["r1", "r2"],
        )
        sparse = ht.mtx_to_sparse_dict(df)
        assert ("r1", "A") in sparse
        assert ("r2", "B") in sparse
        assert ("r1", "B") not in sparse

    def test_prune_identical_reads(self):
        df = pd.DataFrame(
            [[1, 0], [1, 0], [0, 1]],
            columns=["A", "B"],
            index=["r1", "r2", "r3"],
        )
        pruned = ht.prune_identical_reads(df)
        assert pruned.shape[0] == 2

    def test_create_paired_matrix(self):
        b1 = pd.DataFrame(
            [[1, 0], [1, 1]],
            columns=["A", "B"],
            index=["r1", "r2"],
        )
        b2 = pd.DataFrame(
            [[1, 0], [0, 1]],
            columns=["A", "B"],
            index=["r1", "r3"],
        )
        paired, mispaired, unpaired = ht.create_paired_matrix(b1, b2)
        assert "r1" in paired.index  # r1 is in both
        assert "r2" in unpaired.index or "r3" in unpaired.index


class TestPipelineValidation:
    """Tests for pipeline input validation."""

    def test_invalid_file_count(self, tmp_path):
        f1 = tmp_path / "a.fq"
        f2 = tmp_path / "b.fq"
        f3 = tmp_path / "c.fq"
        for f in (f1, f2, f3):
            f.write_text("")
        with pytest.raises(ValueError, match="1.*or 2"):
            run_pipeline([str(f1), str(f2), str(f3)], "dna", str(tmp_path))

    def test_invalid_seq_type(self, tmp_path):
        f = tmp_path / "a.fq"
        f.write_text("")
        with pytest.raises(ValueError, match="dna.*rna"):
            run_pipeline([str(f)], "invalid", str(tmp_path))

    def test_invalid_beta(self, tmp_path):
        f = tmp_path / "a.fq"
        f.write_text("")
        with pytest.raises(ValueError, match="beta"):
            run_pipeline([str(f)], "dna", str(tmp_path), beta=0.5)

    def test_invalid_enumerate(self, tmp_path):
        f = tmp_path / "a.fq"
        f.write_text("")
        with pytest.raises(ValueError, match="enumerate"):
            run_pipeline([str(f)], "dna", str(tmp_path), enumerate_count=0)


@pytest.mark.skipif(
    not (HAS_RAZERS3 and HAS_SOLVER and TEST_RNA_1.exists()),
    reason="Requires razers3, ILP solver, and test data",
)
class TestPipelineIntegration:
    """Integration tests that run the full pipeline (require external tools)."""

    def test_rna_paired_end(self, tmp_path):
        result = run_pipeline(
            input_files=[str(TEST_RNA_1), str(TEST_RNA_2)],
            seq_type="rna",
            output_dir=str(tmp_path),
            verbose=False,
        )
        assert result.result_df is not None
        assert Path(result.output_csv).exists()
        assert Path(result.output_plot).exists()
        # Check result has expected columns
        for col in ["A1", "A2", "B1", "B2", "C1", "C2"]:
            assert col in result.result_4digit.columns
