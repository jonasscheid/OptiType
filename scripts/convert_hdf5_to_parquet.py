#!/usr/bin/env python3
"""
Convert HDF5 allele data to Parquet format.

This script converts the legacy alleles.h5 file to Parquet format
for the modernized OptiType package.

Usage:
    python convert_hdf5_to_parquet.py [--input data/alleles.h5] [--output data/alleles/]
"""

import argparse
import sys
from pathlib import Path

import pandas as pd


def convert_hdf5_to_parquet(input_path: str, output_dir: str) -> None:
    """
    Convert HDF5 allele database to Parquet files.

    Args:
        input_path: Path to input HDF5 file.
        output_dir: Directory for output Parquet files.
    """
    input_path = Path(input_path)
    output_dir = Path(output_dir)

    if not input_path.exists():
        print(f"Error: Input file not found: {input_path}")
        sys.exit(1)

    output_dir.mkdir(parents=True, exist_ok=True)

    print(f"Reading HDF5 file: {input_path}")

    # Read each table using pd.read_hdf which handles the format better
    tables_to_convert = ["table", "features"]

    for table_name in tables_to_convert:
        print(f"Converting: {table_name}")

        try:
            df = pd.read_hdf(str(input_path), table_name)
            print(f"  Shape: {df.shape}")
            print(f"  Columns: {list(df.columns)}")

            output_path = output_dir / f"{table_name}.parquet"
            df.to_parquet(output_path, compression="zstd", index=True)
            print(f"  Saved to: {output_path}")
        except Exception as e:
            print(f"  Error: {e}")

    print("\nConversion complete!")
    print(f"Output directory: {output_dir}")

    # Verify the files
    print("\nVerifying output files:")
    for parquet_file in output_dir.glob("*.parquet"):
        df = pd.read_parquet(parquet_file)
        print(f"  {parquet_file.name}: {df.shape[0]} rows, {df.shape[1]} columns")


def main():
    parser = argparse.ArgumentParser(
        description="Convert OptiType HDF5 allele data to Parquet format"
    )
    parser.add_argument(
        "--input", "-i",
        default="data/alleles.h5",
        help="Path to input HDF5 file (default: data/alleles.h5)",
    )
    parser.add_argument(
        "--output", "-o",
        default="data/alleles/",
        help="Output directory for Parquet files (default: data/alleles/)",
    )

    args = parser.parse_args()
    convert_hdf5_to_parquet(args.input, args.output)


if __name__ == "__main__":
    main()
