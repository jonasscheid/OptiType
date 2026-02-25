[![CI](https://github.com/FRED-2/OptiType/actions/workflows/ci.yml/badge.svg)](https://github.com/FRED-2/OptiType/actions/workflows/ci.yml)
[![PyPI version](https://badge.fury.io/py/optitype.svg)](https://badge.fury.io/py/optitype)

# OptiType

Precision HLA typing from next-generation sequencing data

**Authors:** András Szolek, Benjamin Schubert, Christopher Mohr, Jonas Scheid
**Version:** 1.5.0
**License:** BSD-3-Clause

## Introduction

OptiType is a novel HLA genotyping algorithm based on integer linear
programming, capable of producing accurate 4-digit HLA genotyping predictions
from NGS data by simultaneously selecting all major and minor HLA Class I alleles.

## Installation

### Via pip (Recommended)

```bash
pip install optitype
```

### From source

```bash
git clone https://github.com/FRED-2/OptiType.git
cd OptiType
pip install -e .
```

### External Dependencies

OptiType requires external tools that cannot be installed via pip:

1. **RazerS3** - Read mapper
   ```bash
   # Via conda/bioconda
   conda install -c bioconda razers3

   # Or build from source: https://github.com/seqan/seqan
   ```

2. **ILP Solver** - At least one of:
   - **GLPK** (open source)
     ```bash
     apt install glpk-utils  # Debian/Ubuntu
     conda install -c conda-forge glpk  # Conda
     ```
   - **CBC** (open source)
     ```bash
     apt install coinor-cbc  # Debian/Ubuntu
     conda install -c conda-forge coincbc  # Conda
     ```
   - **CPLEX** (commercial, free for academia)

### Verify Installation

```bash
optitype check-deps
```

## Quick Start

### Command Line

```bash
# DNA sequencing (paired-end)
optitype run -i reads_1.fastq -i reads_2.fastq --dna -o results/

# RNA sequencing (single-end)
optitype run -i sample.fastq --rna -o results/

# Re-analyze from BAM
optitype run -i mapped.bam --dna -o results/
```

### Python API

```python
from optitype import run_hla_typing, HLATypingConfig

result = run_hla_typing(
    fastq_files=["sample_1.fastq", "sample_2.fastq"],
    seq_type="dna",
    config=HLATypingConfig(solver="cbc", threads=4)
)

print(result.best_result)
# {'A1': 'A*02:01', 'A2': 'A*03:01', 'B1': 'B*07:02', ...}
```

## Usage

```
optitype run --help

Usage: optitype run [OPTIONS]

  Run HLA typing analysis.

Options:
  -i, --input PATH       Input FASTQ or BAM files (use multiple times for paired-end)
  -r, --rna              Input data is RNA sequencing
  -d, --dna              Input data is DNA sequencing (default)
  -o, --outdir PATH      Output directory for results (required)
  -p, --prefix TEXT      Output filename prefix (default: timestamp)
  -b, --beta FLOAT       Homozygosity detection parameter (0.0-0.1)
  -e, --enumerate INT    Number of solutions to enumerate
  --solver [glpk|cbc|cplex]  ILP solver to use
  --razers3 PATH         Path to RazerS3 binary
  --threads INT          Number of threads for mapping
  -v, --verbose          Enable verbose output
  -c, --config PATH      Path to config.ini file
  --help                 Show this message and exit
```

### Additional Commands

```bash
# Check dependencies
optitype check-deps

# Generate config file
optitype init-config

# Show installation info
optitype info
```

## Configuration

Generate a config file with:

```bash
optitype init-config -o config.ini
```

Key settings:

```ini
[mapping]
threads=4           # Threads for read mapping

[ilp]
solver=glpk        # ILP solver: glpk, cbc, or cplex
threads=1          # Threads for ILP solver

[behavior]
deletebam=true     # Delete intermediate BAM files
unpaired_weight=0  # Weight for unpaired reads (0-1)
```

## Docker

```bash
docker pull fred2/optitype
docker run -v /path/to/data:/data -t fred2/optitype \
    -i /data/reads_1.fastq -i /data/reads_2.fastq --dna -o /data/results/
```

## Test Examples

```bash
# DNA (paired-end)
optitype run \
    -i ./test/exome/NA11995_SRR766010_1_fished.fastq \
    -i ./test/exome/NA11995_SRR766010_2_fished.fastq \
    --dna -v -o ./test/exome/

# RNA (paired-end)
optitype run \
    -i ./test/rna/CRC_81_N_1_fished.fastq \
    -i ./test/rna/CRC_81_N_2_fished.fastq \
    --rna -v -o ./test/rna/
```

## Output

OptiType produces:
- `*_result.tsv` - HLA typing results
- `*_coverage_plot.pdf` - Coverage visualization

Example output:

```
	A1	A2	B1	B2	C1	C2	Reads	Objective
0	A*02:01	A*03:01	B*07:02	B*44:02	C*07:02	C*05:01	1234	1156.5
```

## Migration from v1.x

Version 1.5 introduces a modernized CLI. Main changes:

- Install with `pip install optitype`
- Use `optitype run` instead of `python OptiTypePipeline.py`
- Multiple input files: use `-i file1 -i file2` instead of `-i file1 file2`
- Data bundled with package (no need to set paths)

The core algorithm and output format remain unchanged.

## Requirements

- Python 3.10+
- External: RazerS3, ILP solver (GLPK/CBC/CPLEX)

## Reference

Szolek, A, Schubert, B, Mohr, C, Sturm, M, Feldhahn, M, and Kohlbacher, O (2014).
**OptiType: precision HLA typing from next-generation sequencing data**
*Bioinformatics*, 30(23):3310-6.
[doi:10.1093/bioinformatics/btu548](https://doi.org/10.1093/bioinformatics/btu548)

## Contact

András Szolek
szolek@informatik.uni-tuebingen.de
University of Tübingen, Applied Bioinformatics
