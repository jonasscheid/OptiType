"""Tests for CLI interface."""

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
    assert "2.0.0" in result.output


def test_cli_run_help():
    """Test that run subcommand help works."""
    runner = CliRunner()
    result = runner.invoke(main, ["run", "--help"])
    assert result.exit_code == 0
    assert "--input" in result.output
    assert "--dna" in result.output
    assert "--rna" in result.output


def test_cli_check_deps():
    """Test that check-deps command works."""
    runner = CliRunner()
    result = runner.invoke(main, ["check-deps"])
    # Exit code may vary based on installed dependencies
    assert "RazerS3" in result.output
    assert "GLPK" in result.output


def test_cli_info():
    """Test that info command works."""
    runner = CliRunner()
    result = runner.invoke(main, ["info"])
    assert result.exit_code == 0
    assert "OptiType version" in result.output
