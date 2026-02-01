"""
Genotype Parser Module
Parses AncestryDNA and 23andMe raw genotype files.
"""

import pandas as pd
import numpy as np
from typing import Dict, List, Tuple, Optional
import re


class GenotypeParser:
    """
    Parser for raw genotype files from AncestryDNA and 23andMe.
    Supports standard formats with RSID, chromosome, position, and genotype data.
    """
    
    def __init__(self):
        self.data = None
        self.metadata = {}
        
    def parse_file(self, filepath: str) -> pd.DataFrame:
        """
        Parse a raw genotype file.
        
        Args:
            filepath: Path to the raw genotype file
            
        Returns:
            DataFrame with columns: rsid, chromosome, position, genotype
        """
        lines = []
        metadata_lines = []
        
        with open(filepath, 'r') as f:
            for line in f:
                line = line.strip()
                if line.startswith('#'):
                    metadata_lines.append(line)
                elif line and not line.startswith('#'):
                    lines.append(line)
        
        # Parse metadata
        self._parse_metadata(metadata_lines)
        
        # Parse data
        if not lines:
            raise ValueError("No data found in file")
        
        # Detect format from first data line
        first_line = lines[0]
        if '\t' in first_line:
            delimiter = '\t'
        else:
            delimiter = ','
        
        # Try to parse as DataFrame
        data_text = '\n'.join(lines)
        from io import StringIO
        
        try:
            self.data = pd.read_csv(StringIO(data_text), sep=delimiter, 
                                   names=['rsid', 'chromosome', 'position', 'genotype'],
                                   dtype={'rsid': str, 'chromosome': str, 
                                          'position': int, 'genotype': str})
        except:
            # Fallback: try with header
            self.data = pd.read_csv(StringIO(data_text), sep=delimiter)
            # Rename columns to standard names
            self.data.columns = ['rsid', 'chromosome', 'position', 'genotype']
        
        # Clean data
        self.data = self._clean_data(self.data)
        
        return self.data
    
    def _parse_metadata(self, metadata_lines: List[str]):
        """Extract metadata from comment lines."""
        for line in metadata_lines:
            if ':' in line:
                key_val = line.lstrip('#').split(':', 1)
                if len(key_val) == 2:
                    self.metadata[key_val[0].strip()] = key_val[1].strip()
    
    def _clean_data(self, df: pd.DataFrame) -> pd.DataFrame:
        """
        Clean and normalize genotype data.
        
        - Remove invalid chromosomes
        - Convert genotypes to standard format
        - Remove no-calls and missing data
        """
        # Keep only autosomal chromosomes and X
        valid_chroms = [str(i) for i in range(1, 23)] + ['X']
        df = df[df['chromosome'].isin(valid_chroms)].copy()
        
        # Remove no-calls and missing data
        df = df[~df['genotype'].isin(['--', '00', 'II', 'DD', 'DI'])].copy()
        df = df.dropna(subset=['genotype'])
        
        # Standardize genotype format (ensure uppercase)
        df['genotype'] = df['genotype'].str.upper()
        
        # Sort by chromosome and position
        df['chrom_num'] = df['chromosome'].replace('X', '23').astype(int)
        df = df.sort_values(['chrom_num', 'position']).drop('chrom_num', axis=1)
        df = df.reset_index(drop=True)
        
        return df
    
    def get_snps(self, rsids: List[str]) -> pd.DataFrame:
        """
        Get specific SNPs by RSID.
        
        Args:
            rsids: List of RSIDs to retrieve
            
        Returns:
            DataFrame with matching SNPs
        """
        if self.data is None:
            raise ValueError("No data loaded. Call parse_file first.")
        
        return self.data[self.data['rsid'].isin(rsids)].copy()
    
    def get_chromosome(self, chromosome: str) -> pd.DataFrame:
        """Get all SNPs from a specific chromosome."""
        if self.data is None:
            raise ValueError("No data loaded. Call parse_file first.")
        
        return self.data[self.data['chromosome'] == str(chromosome)].copy()
    
    def to_numeric_matrix(self) -> Tuple[np.ndarray, List[str]]:
        """
        Convert genotypes to numeric matrix for analysis.
        
        Returns:
            Tuple of (numeric_matrix, rsid_list)
            Matrix encoding: 0=homozygous first, 1=heterozygous, 2=homozygous second
        """
        if self.data is None:
            raise ValueError("No data loaded. Call parse_file first.")
        
        numeric_data = []
        rsids = []
        
        for _, row in self.data.iterrows():
            genotype = row['genotype']
            if len(genotype) == 2:
                allele1, allele2 = genotype[0], genotype[1]
                # Simple encoding based on allele match
                if allele1 == allele2:
                    # Homozygous - use alphabetical order for consistency
                    # A < C < G < T
                    value = 0 if allele1 in ['A', 'C'] else 2
                else:
                    # Heterozygous
                    value = 1
                numeric_data.append(value)
                rsids.append(row['rsid'])
        
        return np.array(numeric_data).reshape(1, -1), rsids
