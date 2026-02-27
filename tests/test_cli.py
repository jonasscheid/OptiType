"""Tests for CLI interface."""

from pathlib import Path
from unittest.mock import patch

from click.testing import CliRunner

from optitype.cli import main


def test_cli_help():
    """Test that help command works."""
    runner = CliRunner()
    result = runner.invoke(main, ["--help"])
    assert result.exit_code == 0
    assert "OptiType" in result.output


def test_cli_version():
    """Test that version command works."""
    runner = CliRunner()
    result = runner.invoke(main, ["--version"])
    assert result.exit_code == 0
    assert "1.5.0" in result.output


def test_cli_run_help():
    """Test that run subcommand help works."""
    runner = CliRunner()
    result = runner.invoke(main, ["run", "--help"])
    assert result.exit_code == 0
    for expected in ["--input", "--dna", "--rna", "--mapper", "--mapper-config", "razers3", "yara"]:
        assert expected in result.output


def test_cli_check_deps():
    """Test that check-deps command works."""
    runner = CliRunner()
    result = runner.invoke(main, ["check-deps"])
    # Exit code may vary based on installed dependencies
    assert "RazerS3" in result.output
    assert "GLPK" in result.output
    assert "YARA" in result.output


def test_cli_info():
    """Test that info command works."""
    runner = CliRunner()
    result = runner.invoke(main, ["info"])
    assert result.exit_code == 0
    assert "OptiType version" in result.output


def test_cli_init_config_mapper_sections(tmp_path):
    """Test that init-config generates mapper config sections."""
    out = tmp_path / "config.ini"
    runner = CliRunner()
    result = runner.invoke(main, ["init-config", "-o", str(out)])
    assert result.exit_code == 0
    content = out.read_text()
    assert "[razers3]" in content
    assert "[yara]" in content
    assert "args =" in content


def test_cli_run_yara_missing_binary(tmp_path):
    """Test that --mapper yara fails gracefully when yara_mapper is not installed."""
    fq = tmp_path / "reads.fq"
    fq.write_text("@read1\nACGT\n+\nIIII\n")
    runner = CliRunner()
    with patch("optitype.cli.shutil.which", side_effect=lambda x: None if x == "yara_mapper" else f"/usr/bin/{x}"):
        result = runner.invoke(
            main,
            ["run", "-i", str(fq), "--dna", "--mapper", "yara", "-o", str(tmp_path / "out")],
        )
    assert result.exit_code != 0
    assert "yara_mapper" in result.output
