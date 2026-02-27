"""
Programmatic API for OptiType.

This module provides a clean Python API for running HLA typing analyses.
"""

from dataclasses import dataclass
from pathlib import Path
from typing import Any

import pandas as pd

from optitype.pipeline import PipelineConfig, load_config, load_mapper_config, run_pipeline


@dataclass
class HLATypingConfig:
    """
    Configuration for HLA typing analysis.

    Attributes:
        solver: ILP solver to use ('glpk', 'cbc', 'cplex').
        threads: Number of threads for ILP solving.
        mapping_threads: Number of threads for read mapping.
        mapper: Read mapper to use ('razers3' or 'yara').
        razers3_path: Path to RazerS3 binary.
        mapper_config: Path to mapper config file with custom alignment args.
        beta: Homozygosity detection parameter (0.0 to 0.1).
        enumerate_count: Number of solutions to enumerate.
        delete_bam: Whether to delete intermediate BAM files.
        unpaired_weight: Weight for unpaired reads in paired-end mode.
        use_discordant: Whether to use discordant read pairs.
    """

    solver: str = "glpk"
    threads: int = 1
    mapping_threads: int = 4
    mapper: str = "razers3"
    razers3_path: str = "razers3"
    mapper_config: str | Path | None = None
    beta: float = 0.009
    enumerate_count: int = 1
    delete_bam: bool = True
    unpaired_weight: float = 0.0
    use_discordant: bool = False

    def to_pipeline_config(self) -> PipelineConfig:
        """Convert to internal PipelineConfig."""
        razers3_args, yara_args = load_mapper_config(self.mapper_config)
        return PipelineConfig(
            razers3_path=self.razers3_path,
            mapper=self.mapper,
            razers3_args=razers3_args,
            yara_args=yara_args,
            mapping_threads=self.mapping_threads,
            solver=self.solver,
            ilp_threads=self.threads,
            delete_bam=self.delete_bam,
            unpaired_weight=self.unpaired_weight,
            use_discordant=self.use_discordant,
        )


@dataclass
class HLATypingResult:
    """
    Results from HLA typing analysis.

    Attributes:
        best_result: Dictionary mapping locus positions to allele types.
        all_results: DataFrame with all enumerated results.
        reads: Number of reads explained by the best solution.
        objective: Objective function value of the best solution.
        output_csv: Path to the output CSV file.
        output_plot: Path to the coverage plot PDF.
    """

    best_result: dict[str, str | None]
    all_results: pd.DataFrame
    reads: float
    objective: float
    output_csv: str | None = None
    output_plot: str | None = None

    def to_dict(self) -> dict[str, Any]:
        """Convert result to a dictionary."""
        return {
            "typing": self.best_result,
            "reads": self.reads,
            "objective": self.objective,
            "output_csv": self.output_csv,
            "output_plot": self.output_plot,
        }

    def __repr__(self) -> str:
        types = [f"{k}: {v}" for k, v in self.best_result.items() if v is not None]
        return f"HLATypingResult({', '.join(types)}, reads={self.reads:.0f})"


def run_hla_typing(
    fastq_files: str | list[str],
    seq_type: str,
    output_dir: str | Path | None = None,
    config: HLATypingConfig | None = None,
    config_file: str | Path | None = None,
    verbose: bool = False,
) -> HLATypingResult:
    """
    Run HLA typing analysis on FASTQ or BAM files.

    This is the main entry point for programmatic use of OptiType.

    Args:
        fastq_files: Path to input file(s). Single file for single-end,
            list of two files for paired-end.
        seq_type: Sequencing type ('dna' or 'rna').
        output_dir: Output directory. If None, uses current directory.
        config: HLATypingConfig with analysis settings. If None, uses defaults.
        config_file: Path to INI config file for additional settings.
        verbose: Enable verbose output.

    Returns:
        HLATypingResult with typing results.

    Raises:
        ValueError: If inputs are invalid.
        FileNotFoundError: If input files or reference data not found.

    Example:
        >>> from optitype import run_hla_typing, HLATypingConfig
        >>> result = run_hla_typing(
        ...     fastq_files=["sample_1.fastq", "sample_2.fastq"],
        ...     seq_type="dna",
        ...     config=HLATypingConfig(solver="cbc", threads=4)
        ... )
        >>> print(result.best_result)
        {'A1': 'A*02:01', 'A2': 'A*03:01', 'B1': 'B*07:02', ...}
    """
    # Normalize input files to list
    if isinstance(fastq_files, str):
        input_files = [fastq_files]
    else:
        input_files = list(fastq_files)

    # Validate inputs
    for f in input_files:
        if not Path(f).exists():
            raise FileNotFoundError(f"Input file not found: {f}")

    # Set up output directory
    if output_dir is None:
        output_dir = Path.cwd() / "optitype_output"
    else:
        output_dir = Path(output_dir)

    # Set up configuration
    if config is None:
        config = HLATypingConfig()

    # Load config file if provided
    pipeline_config = config.to_pipeline_config()
    if config_file:
        file_config = load_config(config_file)
        # Merge file config with API config (API takes precedence)
        if config.razers3_path == "razers3":
            pipeline_config.razers3_path = file_config.razers3_path
        if config.mapper == "razers3":
            pipeline_config.mapper = file_config.mapper
        if config.solver == "glpk":
            pipeline_config.solver = file_config.solver

    # Run pipeline
    result = run_pipeline(
        input_files=input_files,
        seq_type=seq_type,
        output_dir=str(output_dir),
        beta=config.beta,
        enumerate_count=config.enumerate_count,
        config=pipeline_config,
        verbose=verbose,
    )

    # Extract best result
    best_row = result.result_4digit.iloc[0]
    best_result = {
        "A1": best_row.get("A1"),
        "A2": best_row.get("A2"),
        "B1": best_row.get("B1"),
        "B2": best_row.get("B2"),
        "C1": best_row.get("C1"),
        "C2": best_row.get("C2"),
    }

    return HLATypingResult(
        best_result=best_result,
        all_results=result.result_4digit,
        reads=best_row.get("nof_reads", 0),
        objective=best_row.get("obj", 0),
        output_csv=result.output_csv,
        output_plot=result.output_plot,
    )
