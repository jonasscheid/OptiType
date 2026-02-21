"""
Click-based command-line interface for OptiType.

Usage:
    optitype run -i sample_1.fq -i sample_2.fq --dna -o output/
    optitype check-deps
    optitype init-config
"""

import shutil
import sys
from pathlib import Path

import click

from optitype import __version__
from optitype.io.data import get_data_path, get_reference_fasta


@click.group()
@click.version_option(version=__version__, prog_name="optitype")
def main():
    """OptiType: Precision HLA typing from next-generation sequencing data.

    OptiType is a novel HLA genotyping algorithm based on integer linear
    programming, capable of producing accurate 4-digit HLA genotyping
    predictions from NGS data.

    \b
    External dependencies:
      - RazerS3: For read mapping (conda install -c bioconda razers3)
      - ILP solver: GLPK, CBC, or CPLEX

    \b
    Example usage:
      optitype run -i reads_1.fq -i reads_2.fq --dna -o results/
      optitype run -i sample.bam --rna -o results/
    """
    pass


@main.command()
@click.option(
    "-i", "--input",
    "input_files",
    multiple=True,
    required=True,
    type=click.Path(exists=True),
    help="Input FASTQ or BAM files. One file for single-end, two for paired-end.",
)
@click.option(
    "-r", "--rna",
    "seq_type",
    flag_value="rna",
    help="Input data is RNA sequencing.",
)
@click.option(
    "-d", "--dna",
    "seq_type",
    flag_value="dna",
    default=True,
    help="Input data is DNA sequencing (default).",
)
@click.option(
    "-o", "--outdir",
    required=True,
    type=click.Path(),
    help="Output directory for results.",
)
@click.option(
    "-p", "--prefix",
    default=None,
    help="Output filename prefix. Default: timestamp.",
)
@click.option(
    "-b", "--beta",
    type=float,
    default=0.009,
    help="Homozygosity detection parameter (0.0-0.1). Default: 0.009.",
)
@click.option(
    "-e", "--enumerate",
    "enumerate_count",
    type=int,
    default=1,
    help="Number of solutions to enumerate. Default: 1.",
)
@click.option(
    "--solver",
    type=click.Choice(["glpk", "cbc", "cplex"]),
    default="glpk",
    envvar="OPTITYPE_SOLVER",
    help="ILP solver to use. Default: glpk.",
)
@click.option(
    "--razers3",
    type=click.Path(),
    default=None,
    envvar="OPTITYPE_RAZERS3",
    help="Path to RazerS3 binary. Default: search PATH.",
)
@click.option(
    "--threads",
    type=int,
    default=4,
    help="Number of threads for mapping. Default: 4.",
)
@click.option(
    "--ilp-threads",
    type=int,
    default=1,
    help="Number of threads for ILP solver. Default: 1.",
)
@click.option(
    "--keep-bam/--delete-bam",
    default=False,
    help="Keep or delete intermediate BAM files. Default: delete.",
)
@click.option(
    "-c", "--config",
    type=click.Path(exists=True),
    default=None,
    help="Path to config.ini file for additional settings.",
)
@click.option(
    "-v", "--verbose",
    is_flag=True,
    help="Enable verbose output.",
)
def run(
    input_files,
    seq_type,
    outdir,
    prefix,
    beta,
    enumerate_count,
    solver,
    razers3,
    threads,
    ilp_threads,
    keep_bam,
    config,
    verbose,
):
    """Run HLA typing analysis.

    \b
    Examples:
      # Paired-end DNA analysis
      optitype run -i reads_1.fq -i reads_2.fq --dna -o results/

      # Single-end RNA analysis
      optitype run -i sample.fastq --rna -o results/

      # Re-analyze from BAM file
      optitype run -i mapped.bam --dna -o results/

      # With custom settings
      optitype run -i reads.fq --dna -o results/ --solver cbc --threads 8
    """
    from optitype.pipeline import PipelineConfig, run_pipeline

    # Validate inputs
    if len(input_files) not in (1, 2):
        raise click.BadParameter(
            "Number of input files must be 1 (single-end) or 2 (paired-end)",
            param_hint="'-i/--input'",
        )

    if not 0.0 <= beta < 0.1:
        raise click.BadParameter(
            "Beta must be between 0.0 and 0.1",
            param_hint="'-b/--beta'",
        )

    if enumerate_count < 1:
        raise click.BadParameter(
            "Enumerate count must be at least 1",
            param_hint="'-e/--enumerate'",
        )

    # Determine razers3 path
    if razers3 is None:
        razers3 = shutil.which("razers3")
        if razers3 is None:
            raise click.ClickException(
                "RazerS3 not found in PATH. Install with: conda install -c bioconda razers3\n"
                "Or specify path with --razers3 option or OPTITYPE_RAZERS3 environment variable."
            )

    # Create config
    pipeline_config = PipelineConfig(
        razers3_path=razers3,
        mapping_threads=threads,
        solver=solver,
        ilp_threads=ilp_threads,
        delete_bam=not keep_bam,
        unpaired_weight=0.0,
        use_discordant=False,
    )

    # Load additional settings from config file if provided
    if config:
        from optitype.pipeline import load_config
        file_config = load_config(config)
        if file_config.unpaired_weight > 0:
            pipeline_config.unpaired_weight = file_config.unpaired_weight
        if file_config.use_discordant:
            pipeline_config.use_discordant = file_config.use_discordant

    # Run the pipeline
    if verbose:
        click.echo(f"OptiType {__version__}")
        click.echo(f"Input files: {list(input_files)}")
        click.echo(f"Sequence type: {seq_type}")
        click.echo(f"Output directory: {outdir}")
        click.echo(f"Solver: {solver}")
        click.echo()

    try:
        result = run_pipeline(
            input_files=list(input_files),
            seq_type=seq_type,
            output_dir=outdir,
            prefix=prefix,
            beta=beta,
            enumerate_count=enumerate_count,
            config=pipeline_config,
            verbose=verbose,
        )

        # Print results summary
        click.echo()
        click.echo("=" * 60)
        click.echo("HLA Typing Results")
        click.echo("=" * 60)

        best = result.result_4digit.iloc[0]
        for locus in ["A", "B", "C"]:
            a1 = best.get(f"{locus}1", "-")
            a2 = best.get(f"{locus}2", "-")
            if a1 == a2:
                click.echo(f"  HLA-{locus}: {a1} (homozygous)")
            else:
                click.echo(f"  HLA-{locus}: {a1}, {a2}")

        click.echo()
        click.echo(f"Reads: {best['nof_reads']:.0f}")
        click.echo(f"Objective: {best['obj']:.2f}")
        click.echo()
        click.echo(f"Results written to: {result.output_csv}")
        click.echo(f"Coverage plot: {result.output_plot}")

    except Exception as e:
        raise click.ClickException(str(e))


@main.command("check-deps")
def check_deps():
    """Check that all external dependencies are available.

    Verifies that RazerS3 and an ILP solver are installed and accessible.
    """
    click.echo("Checking OptiType dependencies...")
    click.echo()

    all_ok = True

    # Check RazerS3
    razers3 = shutil.which("razers3")
    if razers3:
        click.echo(click.style("  [OK]", fg="green") + f" RazerS3: {razers3}")
    else:
        click.echo(click.style("  [MISSING]", fg="red") + " RazerS3")
        click.echo("    Install with: conda install -c bioconda razers3")
        all_ok = False

    # Check GLPK
    glpsol = shutil.which("glpsol")
    if glpsol:
        click.echo(click.style("  [OK]", fg="green") + f" GLPK: {glpsol}")
    else:
        click.echo(click.style("  [MISSING]", fg="yellow") + " GLPK")
        click.echo("    Install with: apt install glpk-utils  (or conda install -c conda-forge glpk)")

    # Check CBC
    cbc = shutil.which("cbc")
    if cbc:
        click.echo(click.style("  [OK]", fg="green") + f" CBC: {cbc}")
    else:
        click.echo(click.style("  [MISSING]", fg="yellow") + " CBC")
        click.echo("    Install with: apt install coinor-cbc  (or conda install -c conda-forge coincbc)")

    # Check CPLEX (optional)
    cplex = shutil.which("cplex")
    if cplex:
        click.echo(click.style("  [OK]", fg="green") + f" CPLEX: {cplex}")
    else:
        click.echo(click.style("  [N/A]", fg="cyan") + " CPLEX (optional, commercial solver)")

    # Check reference data
    click.echo()
    try:
        data_path = get_data_path()
        click.echo(click.style("  [OK]", fg="green") + f" Reference data: {data_path}")
    except FileNotFoundError:
        click.echo(click.style("  [MISSING]", fg="red") + " Reference data")
        click.echo("    Set OPTITYPE_DATA environment variable to data directory")
        all_ok = False

    click.echo()
    if all_ok and (glpsol or cbc):
        click.echo(click.style("All required dependencies are available!", fg="green"))
    else:
        click.echo(click.style("Some dependencies are missing.", fg="red"))
        if not (glpsol or cbc):
            click.echo("At least one ILP solver (GLPK or CBC) is required.")
        raise SystemExit(1)


@main.command("init-config")
@click.option(
    "-o", "--output",
    type=click.Path(),
    default="config.ini",
    help="Output path for config file. Default: config.ini",
)
@click.option(
    "--force",
    is_flag=True,
    help="Overwrite existing config file.",
)
def init_config(output, force):
    """Generate a default configuration file.

    Creates a config.ini file with default settings that can be customized.
    """
    config_content = """[mapping]

# Path to RazerS3 binary. If not specified, searched in PATH.
# razers3=/path/to/razers3

# Number of threads for read mapping
threads=4

[ilp]

# ILP solver to use: glpk, cbc, or cplex
solver=glpk

# Number of threads for ILP solver
threads=1

[behavior]

# Delete intermediate BAM files after processing
deletebam=true

# Weight for unpaired reads in paired-end mode (0-1)
# 0 means unpaired reads are ignored
# 1 means unpaired reads are weighted equally to paired reads
unpaired_weight=0

# Use discordant read pairs (where ends map to different alleles)
use_discordant=false
"""

    output_path = Path(output)
    if output_path.exists() and not force:
        raise click.ClickException(
            f"Config file already exists: {output}\n"
            "Use --force to overwrite."
        )

    output_path.write_text(config_content)
    click.echo(f"Config file created: {output}")
    click.echo("Edit the file to customize settings.")


@main.command("info")
def info():
    """Show information about the OptiType installation."""
    click.echo(f"OptiType version: {__version__}")
    click.echo()

    # Python info
    click.echo(f"Python: {sys.version}")
    click.echo()

    # Check data path
    try:
        data_path = get_data_path()
        click.echo(f"Data directory: {data_path}")

        # Check for reference files
        for ref_type in ["dna", "rna"]:
            try:
                ref_path = get_reference_fasta(ref_type)
                click.echo(f"  {ref_type.upper()} reference: {ref_path}")
            except FileNotFoundError:
                click.echo(f"  {ref_type.upper()} reference: NOT FOUND")

        # Check for allele data
        alleles_path = data_path / "alleles"
        if alleles_path.exists():
            parquet_files = list(alleles_path.glob("*.parquet"))
            click.echo(f"  Allele data: {len(parquet_files)} parquet files")
        else:
            click.echo("  Allele data: NOT FOUND")
    except FileNotFoundError as e:
        click.echo(f"Data directory: NOT FOUND ({e})")

    click.echo()

    # Package dependencies
    click.echo("Dependencies:")
    for pkg in ["click", "numpy", "pandas", "pyarrow", "pyomo", "pysam", "matplotlib"]:
        try:
            mod = __import__(pkg)
            version = getattr(mod, "__version__", "unknown")
            click.echo(f"  {pkg}: {version}")
        except ImportError:
            click.echo(f"  {pkg}: NOT INSTALLED")


if __name__ == "__main__":
    main()
