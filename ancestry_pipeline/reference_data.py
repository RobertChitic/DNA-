"""
Reference Data Handler Module
Manages European reference genetic datasets for ancestry analysis.
"""

import pandas as pd
import numpy as np
from typing import Dict, List, Optional, Tuple


class ReferenceDataHandler:
    """
    Handles loading and processing of European reference genetic datasets.
    Supports synthetic data generation for testing when real datasets are unavailable.
    """
    
    # European countries and regions for reference
    EUROPEAN_COUNTRIES = [
        'British', 'Irish', 'French', 'Spanish', 'Italian',
        'German', 'Polish', 'Romanian', 'Greek', 'Portuguese',
        'Dutch', 'Belgian', 'Austrian', 'Swiss', 'Swedish',
        'Norwegian', 'Danish', 'Finnish', 'Czech', 'Hungarian',
        'Bulgarian', 'Croatian', 'Serbian', 'Slovenian', 'Slovak',
        'Ukrainian', 'Russian', 'Estonian', 'Latvian', 'Lithuanian'
    ]
    
    def __init__(self):
        self.reference_data = None
        self.reference_labels = None
        self.reference_populations = None
        self.snp_ids = None
        
    def load_reference_data(self, filepath: str, format: str = 'csv') -> Tuple[np.ndarray, List[str]]:
        """
        Load reference genetic data from file.
        
        Args:
            filepath: Path to reference data file
            format: Format of the file ('csv', 'vcf', 'plink')
            
        Returns:
            Tuple of (data_matrix, population_labels)
        """
        if format == 'csv':
            return self._load_csv(filepath)
        elif format == 'vcf':
            raise NotImplementedError("VCF format not yet implemented")
        elif format == 'plink':
            raise NotImplementedError("PLINK format not yet implemented")
        else:
            raise ValueError(f"Unsupported format: {format}")
    
    def _load_csv(self, filepath: str) -> Tuple[np.ndarray, List[str]]:
        """Load reference data from CSV file."""
        df = pd.read_csv(filepath)
        
        # Assume first column is population/country label
        # Remaining columns are SNP values
        self.reference_labels = df.iloc[:, 0].values
        self.reference_data = df.iloc[:, 1:].values
        self.snp_ids = df.columns[1:].tolist()
        
        return self.reference_data, self.reference_labels.tolist()
    
    def generate_synthetic_reference(self, n_samples_per_pop: int = 20, 
                                     n_snps: int = 1000) -> Tuple[np.ndarray, List[str]]:
        """
        Generate synthetic reference data for testing.
        Creates realistic-looking genetic data with population structure.
        
        Args:
            n_samples_per_pop: Number of samples per population
            n_snps: Number of SNPs to generate
            
        Returns:
            Tuple of (data_matrix, population_labels)
        """
        np.random.seed(42)
        
        # Select subset of countries for synthetic data
        countries = self.EUROPEAN_COUNTRIES[:10]  # Use first 10 countries
        
        reference_data = []
        reference_labels = []
        
        # Create population-specific genetic signatures
        # Each population has slightly different allele frequencies
        for i, country in enumerate(countries):
            # Base allele frequencies vary by population
            base_freq = 0.3 + (i / len(countries)) * 0.4  # Range from 0.3 to 0.7
            
            for sample in range(n_samples_per_pop):
                # Generate individual with population-specific tendencies
                individual = np.random.binomial(2, base_freq + np.random.normal(0, 0.1), n_snps)
                individual = np.clip(individual, 0, 2)  # Ensure values are 0, 1, or 2
                
                reference_data.append(individual)
                reference_labels.append(country)
        
        self.reference_data = np.array(reference_data)
        self.reference_labels = np.array(reference_labels)
        self.snp_ids = [f'rs{i}' for i in range(n_snps)]
        
        # Add some population structure (regional clustering)
        self._add_population_structure()
        
        return self.reference_data, self.reference_labels.tolist()
    
    def _add_population_structure(self):
        """Add realistic population structure to synthetic data."""
        # Define regional clusters
        regions = {
            'Western': ['British', 'Irish', 'French', 'Dutch', 'Belgian'],
            'Southern': ['Spanish', 'Italian', 'Portuguese', 'Greek'],
            'Northern': ['Swedish', 'Norwegian', 'Danish', 'Finnish'],
            'Eastern': ['Polish', 'Romanian', 'Czech', 'Hungarian']
        }
        
        # Adjust data to create regional structure
        for region, countries in regions.items():
            region_indices = [i for i, label in enumerate(self.reference_labels) 
                            if label in countries]
            if region_indices:
                # Add small regional signal to certain SNPs
                region_snps = np.random.choice(len(self.snp_ids), 
                                             size=len(self.snp_ids) // 10, 
                                             replace=False)
                self.reference_data[region_indices][:, region_snps] += np.random.choice(
                    [-1, 1], size=(len(region_indices), len(region_snps))
                )
                self.reference_data = np.clip(self.reference_data, 0, 2)
    
    def get_population_centroids(self) -> Dict[str, np.ndarray]:
        """
        Calculate centroid (mean) for each population.
        
        Returns:
            Dictionary mapping population names to centroid vectors
        """
        if self.reference_data is None or self.reference_labels is None:
            raise ValueError("No reference data loaded")
        
        centroids = {}
        unique_pops = np.unique(self.reference_labels)
        
        for pop in unique_pops:
            pop_indices = self.reference_labels == pop
            centroids[pop] = np.mean(self.reference_data[pop_indices], axis=0)
        
        return centroids
    
    def align_snps(self, user_snps: List[str]) -> Tuple[List[int], List[int]]:
        """
        Align user SNPs with reference SNPs.
        
        Args:
            user_snps: List of SNP IDs from user data
            
        Returns:
            Tuple of (user_indices, reference_indices) for matching SNPs
        """
        if self.snp_ids is None:
            raise ValueError("No reference data loaded")
        
        # Find intersection of SNPs
        reference_set = set(self.snp_ids)
        user_indices = []
        reference_indices = []
        
        for i, snp in enumerate(user_snps):
            if snp in reference_set:
                user_indices.append(i)
                reference_indices.append(self.snp_ids.index(snp))
        
        return user_indices, reference_indices
    
    def get_reference_stats(self) -> Dict:
        """Get statistics about the reference dataset."""
        if self.reference_data is None:
            return {}
        
        unique_pops, counts = np.unique(self.reference_labels, return_counts=True)
        
        return {
            'n_samples': len(self.reference_data),
            'n_snps': self.reference_data.shape[1] if len(self.reference_data.shape) > 1 else 0,
            'n_populations': len(unique_pops),
            'populations': dict(zip(unique_pops, counts.tolist()))
        }
