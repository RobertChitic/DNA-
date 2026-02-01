#!/usr/bin/env python3
"""
Script: 02_convert_ancestry_file.py
Description: Convert AncestryDNA raw genotype file to PLINK binary format

AncestryDNA format:
    rsid    chromosome    position    allele1    allele2
    rs1234  1            12345       A          G

Output: PLINK binary files (.bed, .bim, .fam)

Usage:
    python 02_convert_ancestry_file.py --input <ancestry_file.txt> --output <output_prefix>
"""

import argparse
import os
import sys
import subprocess
import tempfile
import logging
from pathlib import Path
from typing import Optional, Tuple

import pandas as pd
import numpy as np

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)


def validate_input_file(filepath: str) -> bool:
    """Validate that the input file exists and has expected format."""
    if not os.path.exists(filepath):
        logger.error(f"Input file not found: {filepath}")
        return False
    
    # Check file format
    try:
        with open(filepath, 'r') as f:
            # Skip comment lines (AncestryDNA files start with #)
            header = None
            for line in f:
                if not line.startswith('#'):
                    header = line.strip().lower()
                    break
            
            if header is None:
                logger.error("File appears to be empty or contains only comments")
                return False
            
            # Check for expected columns
            expected_cols = {'rsid', 'chromosome', 'position', 'allele1', 'allele2'}
            actual_cols = set(header.split('\t'))
            
            if not expected_cols.issubset(actual_cols):
                # Also check for space-delimited format
                actual_cols = set(header.split())
                if not expected_cols.issubset(actual_cols):
                    logger.warning(f"Non-standard column names detected: {actual_cols}")
                    logger.info("Will attempt to parse file anyway...")
                    
    except Exception as e:
        logger.error(f"Error reading input file: {e}")
        return False
    
    return True


def read_ancestry_file(filepath: str) -> pd.DataFrame:
    """
    Read AncestryDNA raw genotype file.
    
    Handles various formats and comment lines.
    """
    logger.info(f"Reading AncestryDNA file: {filepath}")
    
    # First, determine delimiter and skip rows
    skip_rows = 0
    delimiter = '\t'
    
    with open(filepath, 'r') as f:
        for i, line in enumerate(f):
            if line.startswith('#'):
                skip_rows = i + 1
            else:
                # Detect delimiter
                if '\t' in line:
                    delimiter = '\t'
                else:
                    delimiter = r'\s+'
                break
    
    # Read the file
    df = pd.read_csv(
        filepath,
        sep=delimiter,
        skiprows=skip_rows,
        dtype={
            'rsid': str,
            'chromosome': str,
            'position': int,
            'allele1': str,
            'allele2': str
        },
        comment='#',
        na_values=['--', '-', '0', 'N/A', 'NA', '.']
    )
    
    # Normalize column names
    df.columns = df.columns.str.lower().str.strip()
    
    # Rename columns if needed
    col_mapping = {
        'rsid': 'rsid',
        'rs': 'rsid',
        'snp': 'rsid',
        'marker': 'rsid',
        'chr': 'chromosome',
        'chrom': 'chromosome',
        'chromosome': 'chromosome',
        'pos': 'position',
        'position': 'position',
        'physical_position': 'position',
        'a1': 'allele1',
        'allele1': 'allele1',
        'genotype1': 'allele1',
        'a2': 'allele2',
        'allele2': 'allele2',
        'genotype2': 'allele2'
    }
    
    df = df.rename(columns={col: col_mapping.get(col, col) for col in df.columns})
    
    logger.info(f"Read {len(df):,} variants from input file")
    
    return df


def clean_genotype_data(df: pd.DataFrame) -> pd.DataFrame:
    """
    Clean and standardize genotype data.
    
    - Remove invalid chromosomes
    - Handle missing data
    - Standardize allele codes
    """
    logger.info("Cleaning genotype data...")
    
    original_count = len(df)
    
    # Remove rows with missing essential data
    df = df.dropna(subset=['rsid', 'chromosome', 'position'])
    logger.info(f"Removed {original_count - len(df):,} rows with missing essential data")
    
    # Filter valid chromosomes (1-22, X, Y, MT)
    valid_chroms = [str(i) for i in range(1, 23)] + ['X', 'Y', 'MT', 'M']
    df['chromosome'] = df['chromosome'].astype(str).str.upper().str.replace('CHR', '')
    df = df[df['chromosome'].isin(valid_chroms)]
    logger.info(f"After chromosome filter: {len(df):,} variants")
    
    # Clean alleles (uppercase, handle no-calls)
    df['allele1'] = df['allele1'].fillna('0').astype(str).str.upper()
    df['allele2'] = df['allele2'].fillna('0').astype(str).str.upper()
    
    # Replace invalid alleles with missing
    valid_alleles = {'A', 'C', 'G', 'T', '0', '-', 'I', 'D'}
    df.loc[~df['allele1'].isin(valid_alleles), 'allele1'] = '0'
    df.loc[~df['allele2'].isin(valid_alleles), 'allele2'] = '0'
    
    # Convert position to integer
    df['position'] = df['position'].astype(int)
    
    # Sort by chromosome and position
    chrom_order = {str(i): i for i in range(1, 23)}
    chrom_order.update({'X': 23, 'Y': 24, 'MT': 25, 'M': 25})
    df['chrom_num'] = df['chromosome'].map(chrom_order)
    df = df.sort_values(['chrom_num', 'position']).drop(columns=['chrom_num'])
    
    # Remove duplicates (keep first occurrence)
    df = df.drop_duplicates(subset=['rsid'], keep='first')
    
    logger.info(f"Final variant count: {len(df):,}")
    
    return df


def write_plink_files(df: pd.DataFrame, output_prefix: str, sample_id: str = "USER001") -> bool:
    """
    Write PLINK format files.
    
    Creates:
    - .ped file (genotype data)
    - .map file (variant information)
    
    Then converts to binary format using PLINK.
    """
    logger.info(f"Writing PLINK files with prefix: {output_prefix}")
    
    # Ensure output directory exists
    output_dir = os.path.dirname(output_prefix)
    if output_dir:
        os.makedirs(output_dir, exist_ok=True)
    
    # Write MAP file
    map_file = f"{output_prefix}.map"
    logger.info(f"Writing MAP file: {map_file}")
    
    # PLINK MAP format: chromosome, variant_id, genetic_distance, position
    map_df = df[['chromosome', 'rsid', 'position']].copy()
    map_df['genetic_distance'] = 0  # Placeholder
    map_df = map_df[['chromosome', 'rsid', 'genetic_distance', 'position']]
    
    # Convert chromosome codes for PLINK
    # PLINK uses numeric codes: 23=X, 24=Y, 26=MT (mitochondrial)
    chrom_map = {'X': '23', 'Y': '24', 'MT': '26', 'M': '26'}
    map_df['chromosome'] = map_df['chromosome'].replace(chrom_map)
    
    map_df.to_csv(map_file, sep='\t', header=False, index=False)
    
    # Write PED file
    ped_file = f"{output_prefix}.ped"
    logger.info(f"Writing PED file: {ped_file}")
    
    # PED format: FID IID PID MID SEX PHENOTYPE GENOTYPES...
    # Genotypes are space-separated pairs of alleles
    with open(ped_file, 'w') as f:
        # Family ID, Individual ID, Paternal ID, Maternal ID, Sex (0=unknown), Phenotype (-9=missing)
        f.write(f"{sample_id}\t{sample_id}\t0\t0\t0\t-9")
        
        # Write genotypes
        for _, row in df.iterrows():
            a1 = row['allele1'] if row['allele1'] not in ['0', '-', 'N'] else '0'
            a2 = row['allele2'] if row['allele2'] not in ['0', '-', 'N'] else '0'
            f.write(f"\t{a1}\t{a2}")
        
        f.write("\n")
    
    logger.info("Converting to PLINK binary format...")
    
    # Convert to binary format using PLINK
    try:
        cmd = [
            'plink',
            '--file', output_prefix,
            '--make-bed',
            '--out', output_prefix,
            '--allow-no-sex',
            '--allow-extra-chr'
        ]
        
        result = subprocess.run(
            cmd,
            capture_output=True,
            text=True,
            check=True
        )
        
        logger.info("Successfully converted to PLINK binary format")
        
        # Clean up temporary text files (optional)
        # os.remove(ped_file)
        # os.remove(map_file)
        
        return True
        
    except subprocess.CalledProcessError as e:
        logger.error(f"PLINK conversion failed: {e.stderr}")
        return False
    except FileNotFoundError:
        logger.error("PLINK not found. Please ensure PLINK is installed and in PATH.")
        return False


def get_summary_stats(df: pd.DataFrame) -> dict:
    """Generate summary statistics for the genotype data."""
    stats = {
        'total_variants': len(df),
        'chromosomes': df['chromosome'].nunique(),
        'variant_breakdown': {}
    }
    
    # Count by chromosome
    stats['variant_breakdown'] = df['chromosome'].value_counts().to_dict()
    
    # Missing rate
    missing_a1 = (df['allele1'] == '0').sum()
    missing_a2 = (df['allele2'] == '0').sum()
    stats['missing_rate'] = (missing_a1 + missing_a2) / (len(df) * 2) * 100
    
    return stats


def main():
    parser = argparse.ArgumentParser(
        description='Convert AncestryDNA raw genotype file to PLINK binary format',
        formatter_class=argparse.RawDescriptionHelpFormatter
    )
    
    parser.add_argument(
        '--input', '-i',
        required=True,
        help='Path to AncestryDNA raw data file'
    )
    
    parser.add_argument(
        '--output', '-o',
        required=True,
        help='Output prefix for PLINK files (e.g., data/processed/user)'
    )
    
    parser.add_argument(
        '--sample-id', '-s',
        default='USER001',
        help='Sample ID to use in PLINK files (default: USER001)'
    )
    
    parser.add_argument(
        '--verbose', '-v',
        action='store_true',
        help='Enable verbose output'
    )
    
    args = parser.parse_args()
    
    if args.verbose:
        logger.setLevel(logging.DEBUG)
    
    # Validate input
    if not validate_input_file(args.input):
        sys.exit(1)
    
    # Read and process data
    try:
        df = read_ancestry_file(args.input)
        df = clean_genotype_data(df)
        
        # Print summary
        stats = get_summary_stats(df)
        logger.info(f"Summary Statistics:")
        logger.info(f"  Total variants: {stats['total_variants']:,}")
        logger.info(f"  Chromosomes: {stats['chromosomes']}")
        logger.info(f"  Missing rate: {stats['missing_rate']:.2f}%")
        
        # Write PLINK files
        success = write_plink_files(df, args.output, args.sample_id)
        
        if success:
            logger.info("Conversion completed successfully!")
            logger.info(f"Output files:")
            logger.info(f"  {args.output}.bed")
            logger.info(f"  {args.output}.bim")
            logger.info(f"  {args.output}.fam")
            sys.exit(0)
        else:
            logger.error("Conversion failed!")
            sys.exit(1)
            
    except Exception as e:
        logger.error(f"Error processing file: {e}")
        sys.exit(1)


if __name__ == '__main__':
    main()
