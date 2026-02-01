#!/usr/bin/env python3
"""
Script: 06_pca_probabilities.py
Description: Calculate ancestry probabilities from PCA results using Gaussian kernel

Method:
1. Load PCA coordinates (eigenvec file)
2. Calculate Euclidean distances between user and reference samples
3. Convert distances to probabilities using Gaussian kernel
4. Aggregate by country/population

Usage:
    python 06_pca_probabilities.py --eigenvec <file> --metadata <file> --output <file>
"""

import argparse
import os
import sys
import logging
from typing import Dict, List, Tuple, Optional

import numpy as np
import pandas as pd
from scipy.spatial.distance import cdist
from collections import defaultdict

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)

# Default parameters
DEFAULT_N_PCS = 10
DEFAULT_BANDWIDTH = 0.5


def load_eigenvec(filepath: str, n_pcs: int = DEFAULT_N_PCS) -> pd.DataFrame:
    """
    Load PLINK eigenvec file.
    
    Format: FID IID PC1 PC2 ... PCn
    """
    logger.info(f"Loading eigenvec file: {filepath}")
    
    # Try to detect if file has header
    with open(filepath, 'r') as f:
        first_line = f.readline().strip()
        has_header = 'PC1' in first_line.upper() or 'FID' in first_line.upper()
    
    if has_header:
        df = pd.read_csv(filepath, sep=r'\s+', header=0)
    else:
        # Generate column names
        cols = ['FID', 'IID'] + [f'PC{i}' for i in range(1, 100)]
        df = pd.read_csv(filepath, sep=r'\s+', header=None, names=cols)
    
    # Normalize column names
    df.columns = df.columns.str.upper()
    
    # Rename first two columns if needed
    if df.columns[0] != 'FID':
        df = df.rename(columns={df.columns[0]: 'FID', df.columns[1]: 'IID'})
    
    # Select only required PCs
    pc_cols = [f'PC{i}' for i in range(1, n_pcs + 1)]
    available_pcs = [col for col in pc_cols if col in df.columns]
    
    if len(available_pcs) < n_pcs:
        logger.warning(f"Only {len(available_pcs)} PCs available, requested {n_pcs}")
    
    df = df[['FID', 'IID'] + available_pcs]
    
    logger.info(f"Loaded {len(df)} samples with {len(available_pcs)} PCs")
    
    return df


def load_metadata(filepath: str) -> pd.DataFrame:
    """
    Load sample metadata file.
    
    Expected columns: sample_id, population, country
    """
    logger.info(f"Loading metadata file: {filepath}")
    
    df = pd.read_csv(filepath, sep='\t')
    df.columns = df.columns.str.lower()
    
    # Rename columns if needed
    if 'sample_id' not in df.columns and 'iid' in df.columns:
        df = df.rename(columns={'iid': 'sample_id'})
    
    required_cols = ['sample_id', 'population', 'country']
    for col in required_cols:
        if col not in df.columns:
            logger.warning(f"Missing column '{col}' in metadata")
    
    logger.info(f"Loaded metadata for {len(df)} samples")
    
    return df


def identify_user_samples(eigenvec_df: pd.DataFrame, metadata_df: pd.DataFrame) -> List[str]:
    """
    Identify user samples (not in reference metadata).
    """
    if metadata_df is None or metadata_df.empty:
        # Assume first sample is user
        return [eigenvec_df['IID'].iloc[0]]
    
    reference_ids = set(metadata_df['sample_id'].astype(str))
    user_samples = []
    
    for iid in eigenvec_df['IID']:
        if str(iid) not in reference_ids:
            user_samples.append(iid)
    
    logger.info(f"Identified {len(user_samples)} user sample(s)")
    
    return user_samples


def gaussian_kernel(distances: np.ndarray, bandwidth: float) -> np.ndarray:
    """
    Convert distances to similarities using Gaussian kernel.
    
    similarity = exp(-distance^2 / (2 * bandwidth^2))
    """
    return np.exp(-distances**2 / (2 * bandwidth**2))


def calculate_ancestry_probabilities(
    user_coords: np.ndarray,
    ref_coords: np.ndarray,
    ref_populations: List[str],
    ref_countries: List[str],
    bandwidth: float = DEFAULT_BANDWIDTH
) -> Tuple[Dict[str, float], Dict[str, float]]:
    """
    Calculate ancestry probabilities using Gaussian kernel method.
    
    Returns:
        - Population-level probabilities
        - Country-level probabilities
    """
    # Calculate Euclidean distances
    distances = cdist(user_coords.reshape(1, -1), ref_coords, metric='euclidean')[0]
    
    # Convert to similarities
    similarities = gaussian_kernel(distances, bandwidth)
    
    # Normalize
    similarities = similarities / similarities.sum()
    
    # Aggregate by population
    pop_probs = defaultdict(float)
    for pop, sim in zip(ref_populations, similarities):
        pop_probs[pop] += sim
    
    # Aggregate by country
    country_probs = defaultdict(float)
    for country, sim in zip(ref_countries, similarities):
        country_probs[country] += sim
    
    return dict(pop_probs), dict(country_probs)


def bootstrap_confidence_intervals(
    user_coords: np.ndarray,
    ref_coords: np.ndarray,
    ref_populations: List[str],
    ref_countries: List[str],
    bandwidth: float = DEFAULT_BANDWIDTH,
    n_bootstrap: int = 1000,
    confidence_level: float = 0.95
) -> Dict[str, Tuple[float, float]]:
    """
    Calculate confidence intervals using bootstrap resampling.
    """
    logger.info(f"Running {n_bootstrap} bootstrap iterations...")
    
    n_ref = len(ref_coords)
    country_samples = defaultdict(list)
    
    for i in range(n_bootstrap):
        # Resample reference with replacement
        indices = np.random.choice(n_ref, size=n_ref, replace=True)
        
        sampled_coords = ref_coords[indices]
        sampled_countries = [ref_countries[j] for j in indices]
        
        # Calculate distances and probabilities
        distances = cdist(user_coords.reshape(1, -1), sampled_coords, metric='euclidean')[0]
        similarities = gaussian_kernel(distances, bandwidth)
        similarities = similarities / similarities.sum()
        
        # Aggregate
        country_probs = defaultdict(float)
        for country, sim in zip(sampled_countries, similarities):
            country_probs[country] += sim
        
        for country, prob in country_probs.items():
            country_samples[country].append(prob)
    
    # Calculate confidence intervals
    alpha = (1 - confidence_level) / 2
    ci = {}
    
    for country, samples in country_samples.items():
        lower = np.percentile(samples, alpha * 100)
        upper = np.percentile(samples, (1 - alpha) * 100)
        ci[country] = (lower, upper)
    
    return ci


def main():
    parser = argparse.ArgumentParser(
        description='Calculate ancestry probabilities from PCA results',
        formatter_class=argparse.RawDescriptionHelpFormatter
    )
    
    parser.add_argument(
        '--eigenvec', '-e',
        required=True,
        help='Path to PLINK eigenvec file'
    )
    
    parser.add_argument(
        '--metadata', '-m',
        required=True,
        help='Path to sample metadata file'
    )
    
    parser.add_argument(
        '--output', '-o',
        required=True,
        help='Output file for ancestry probabilities'
    )
    
    parser.add_argument(
        '--n-pcs', '-n',
        type=int,
        default=DEFAULT_N_PCS,
        help=f'Number of PCs to use (default: {DEFAULT_N_PCS})'
    )
    
    parser.add_argument(
        '--bandwidth', '-b',
        type=float,
        default=DEFAULT_BANDWIDTH,
        help=f'Gaussian kernel bandwidth (default: {DEFAULT_BANDWIDTH})'
    )
    
    parser.add_argument(
        '--bootstrap', '-B',
        type=int,
        default=0,
        help='Number of bootstrap iterations for CI (0 to disable)'
    )
    
    parser.add_argument(
        '--user-id',
        help='Specific user sample ID (default: auto-detect)'
    )
    
    args = parser.parse_args()
    
    # Load data
    eigenvec_df = load_eigenvec(args.eigenvec, args.n_pcs)
    metadata_df = load_metadata(args.metadata)
    
    # Identify user sample
    if args.user_id:
        user_samples = [args.user_id]
    else:
        user_samples = identify_user_samples(eigenvec_df, metadata_df)
    
    if not user_samples:
        logger.error("No user samples found!")
        sys.exit(1)
    
    # Merge eigenvec with metadata
    merged = eigenvec_df.merge(
        metadata_df, 
        left_on='IID', 
        right_on='sample_id', 
        how='left'
    )
    
    # Separate user and reference samples
    pc_cols = [col for col in merged.columns if col.startswith('PC')]
    
    user_mask = merged['IID'].isin(user_samples)
    ref_mask = ~user_mask & merged['population'].notna()
    
    user_data = merged[user_mask][pc_cols].values
    ref_data = merged[ref_mask][pc_cols].values
    ref_pops = merged[ref_mask]['population'].tolist()
    ref_countries = merged[ref_mask]['country'].tolist()
    
    logger.info(f"User samples: {len(user_data)}")
    logger.info(f"Reference samples: {len(ref_data)}")
    
    # Calculate probabilities for each user
    results = []
    
    for i, user_id in enumerate(user_samples):
        logger.info(f"Processing user: {user_id}")
        
        pop_probs, country_probs = calculate_ancestry_probabilities(
            user_data[i],
            ref_data,
            ref_pops,
            ref_countries,
            args.bandwidth
        )
        
        # Bootstrap CI if requested
        ci = {}
        if args.bootstrap > 0:
            ci = bootstrap_confidence_intervals(
                user_data[i],
                ref_data,
                ref_pops,
                ref_countries,
                args.bandwidth,
                args.bootstrap
            )
        
        # Store results
        for country, prob in sorted(country_probs.items(), key=lambda x: -x[1]):
            result = {
                'User': user_id,
                'Country': country,
                'PCA_Probability': prob * 100
            }
            
            if ci and country in ci:
                result['CI_Lower'] = ci[country][0] * 100
                result['CI_Upper'] = ci[country][1] * 100
            
            results.append(result)
    
    # Save results
    results_df = pd.DataFrame(results)
    
    output_dir = os.path.dirname(args.output)
    if output_dir:
        os.makedirs(output_dir, exist_ok=True)
    
    results_df.to_csv(args.output, index=False)
    logger.info(f"Results saved to: {args.output}")
    
    # Print summary
    print("\n=== PCA-Based Ancestry Probabilities ===\n")
    for user_id in user_samples:
        user_results = results_df[results_df['User'] == user_id]
        print(f"User: {user_id}")
        print("-" * 50)
        for _, row in user_results.head(10).iterrows():
            prob = row['PCA_Probability']
            country = row['Country']
            if 'CI_Lower' in row:
                print(f"  {country}: {prob:.1f}% ({row['CI_Lower']:.1f}% - {row['CI_Upper']:.1f}%)")
            else:
                print(f"  {country}: {prob:.1f}%")
        print()


if __name__ == '__main__':
    main()
