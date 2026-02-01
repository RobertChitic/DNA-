#!/usr/bin/env python3
"""
Script: 07_admixture_probabilities.py
Description: Calculate ancestry probabilities from ADMIXTURE results

Method:
1. Load ADMIXTURE Q file (individual ancestry proportions)
2. Calculate mean component profiles for each reference country
3. Calculate cosine similarity between user and reference profiles
4. Convert similarities to probabilities

Usage:
    python 07_admixture_probabilities.py --q-file <file> --metadata <file> --output <file>
"""

import argparse
import os
import sys
import logging
from typing import Dict, List, Tuple, Optional

import numpy as np
import pandas as pd
from scipy.spatial.distance import cosine
from collections import defaultdict

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)


def load_q_file(filepath: str, fam_filepath: str) -> pd.DataFrame:
    """
    Load ADMIXTURE Q file.
    
    Q file format: One row per individual, one column per ancestral component
    No sample IDs - must be joined with FAM file.
    """
    logger.info(f"Loading Q file: {filepath}")
    
    # Load Q file (no header)
    q_data = pd.read_csv(filepath, sep=r'\s+', header=None)
    n_components = q_data.shape[1]
    q_data.columns = [f'K{i+1}' for i in range(n_components)]
    
    logger.info(f"Loaded {len(q_data)} samples with {n_components} components")
    
    # Load FAM file for sample IDs
    logger.info(f"Loading FAM file: {fam_filepath}")
    fam_data = pd.read_csv(fam_filepath, sep=r'\s+', header=None, 
                          names=['FID', 'IID', 'PID', 'MID', 'Sex', 'Phenotype'])
    
    if len(fam_data) != len(q_data):
        logger.warning(f"FAM ({len(fam_data)}) and Q ({len(q_data)}) row count mismatch!")
        # Trim to minimum
        min_len = min(len(fam_data), len(q_data))
        fam_data = fam_data.iloc[:min_len]
        q_data = q_data.iloc[:min_len]
    
    # Combine
    result = pd.concat([fam_data[['FID', 'IID']], q_data], axis=1)
    
    return result


def load_metadata(filepath: str) -> pd.DataFrame:
    """
    Load sample metadata file.
    """
    logger.info(f"Loading metadata file: {filepath}")
    
    df = pd.read_csv(filepath, sep='\t')
    df.columns = df.columns.str.lower()
    
    if 'sample_id' not in df.columns and 'iid' in df.columns:
        df = df.rename(columns={'iid': 'sample_id'})
    
    logger.info(f"Loaded metadata for {len(df)} samples")
    
    return df


def identify_user_samples(q_df: pd.DataFrame, metadata_df: pd.DataFrame) -> List[str]:
    """
    Identify user samples (not in reference metadata).
    """
    if metadata_df is None or metadata_df.empty:
        return [q_df['IID'].iloc[0]]
    
    reference_ids = set(metadata_df['sample_id'].astype(str))
    user_samples = []
    
    for iid in q_df['IID']:
        if str(iid) not in reference_ids:
            user_samples.append(iid)
    
    logger.info(f"Identified {len(user_samples)} user sample(s)")
    
    return user_samples


def calculate_country_profiles(
    q_df: pd.DataFrame,
    metadata_df: pd.DataFrame,
    k_cols: List[str]
) -> Dict[str, np.ndarray]:
    """
    Calculate mean ADMIXTURE profile for each country.
    """
    # Merge with metadata
    merged = q_df.merge(
        metadata_df,
        left_on='IID',
        right_on='sample_id',
        how='left'
    )
    
    # Calculate mean profile per country
    profiles = {}
    
    for country in merged['country'].dropna().unique():
        country_data = merged[merged['country'] == country][k_cols].values
        if len(country_data) > 0:
            profiles[country] = np.mean(country_data, axis=0)
            logger.debug(f"Country '{country}': {len(country_data)} samples")
    
    logger.info(f"Calculated profiles for {len(profiles)} countries")
    
    return profiles


def cosine_similarity(a: np.ndarray, b: np.ndarray) -> float:
    """
    Calculate cosine similarity between two vectors.
    """
    return 1 - cosine(a, b)


def calculate_ancestry_probabilities(
    user_profile: np.ndarray,
    country_profiles: Dict[str, np.ndarray]
) -> Dict[str, float]:
    """
    Calculate ancestry probabilities using cosine similarity.
    """
    similarities = {}
    
    for country, profile in country_profiles.items():
        sim = cosine_similarity(user_profile, profile)
        # Ensure non-negative
        similarities[country] = max(0, sim)
    
    # Normalize to probabilities
    total = sum(similarities.values())
    if total > 0:
        probabilities = {k: v / total for k, v in similarities.items()}
    else:
        # Fall back to uniform
        n = len(similarities)
        probabilities = {k: 1/n for k in similarities}
    
    return probabilities


def interpret_components(
    country_profiles: Dict[str, np.ndarray],
    k: int
) -> Dict[int, List[str]]:
    """
    Identify which countries load highest on each component.
    """
    interpretations = {}
    
    for i in range(k):
        # Get component values for each country
        country_vals = [(country, profile[i]) for country, profile in country_profiles.items()]
        # Sort by value
        sorted_countries = sorted(country_vals, key=lambda x: -x[1])
        # Keep top contributors
        top_countries = [c for c, v in sorted_countries if v > 0.1][:3]
        interpretations[i + 1] = top_countries if top_countries else ['Mixed']
    
    return interpretations


def bootstrap_confidence_intervals(
    user_profile: np.ndarray,
    ref_profiles: np.ndarray,
    ref_countries: List[str],
    n_bootstrap: int = 1000,
    confidence_level: float = 0.95
) -> Dict[str, Tuple[float, float]]:
    """
    Calculate confidence intervals using bootstrap resampling of reference samples.
    """
    logger.info(f"Running {n_bootstrap} bootstrap iterations...")
    
    country_samples = defaultdict(list)
    unique_countries = list(set(ref_countries))
    
    for _ in range(n_bootstrap):
        # Bootstrap country profiles
        boot_profiles = {}
        
        for country in unique_countries:
            country_mask = np.array([c == country for c in ref_countries])
            country_data = ref_profiles[country_mask]
            
            if len(country_data) > 0:
                # Sample with replacement
                indices = np.random.choice(len(country_data), size=len(country_data), replace=True)
                boot_profiles[country] = np.mean(country_data[indices], axis=0)
        
        # Calculate probabilities
        probs = calculate_ancestry_probabilities(user_profile, boot_profiles)
        
        for country, prob in probs.items():
            country_samples[country].append(prob)
    
    # Calculate CIs
    alpha = (1 - confidence_level) / 2
    ci = {}
    
    for country, samples in country_samples.items():
        ci[country] = (
            np.percentile(samples, alpha * 100),
            np.percentile(samples, (1 - alpha) * 100)
        )
    
    return ci


def main():
    parser = argparse.ArgumentParser(
        description='Calculate ancestry probabilities from ADMIXTURE results',
        formatter_class=argparse.RawDescriptionHelpFormatter
    )
    
    parser.add_argument(
        '--q-file', '-q',
        required=True,
        help='Path to ADMIXTURE Q file'
    )
    
    parser.add_argument(
        '--fam-file', '-f',
        required=True,
        help='Path to PLINK FAM file (for sample IDs)'
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
    q_df = load_q_file(args.q_file, args.fam_file)
    metadata_df = load_metadata(args.metadata)
    
    # Get component columns
    k_cols = [col for col in q_df.columns if col.startswith('K')]
    k = len(k_cols)
    logger.info(f"Number of ancestral components (K): {k}")
    
    # Identify user samples
    if args.user_id:
        user_samples = [args.user_id]
    else:
        user_samples = identify_user_samples(q_df, metadata_df)
    
    if not user_samples:
        logger.error("No user samples found!")
        sys.exit(1)
    
    # Merge data
    merged = q_df.merge(
        metadata_df,
        left_on='IID',
        right_on='sample_id',
        how='left'
    )
    
    # Calculate country profiles
    country_profiles = calculate_country_profiles(q_df, metadata_df, k_cols)
    
    # Get component interpretations
    component_interp = interpret_components(country_profiles, k)
    
    # Prepare for bootstrap CI
    ref_mask = merged['country'].notna()
    ref_profiles = merged[ref_mask][k_cols].values
    ref_countries = merged[ref_mask]['country'].tolist()
    
    # Calculate probabilities for each user
    results = []
    user_components = {}
    
    for user_id in user_samples:
        logger.info(f"Processing user: {user_id}")
        
        user_data = q_df[q_df['IID'] == user_id][k_cols].values
        
        if len(user_data) == 0:
            logger.warning(f"User {user_id} not found in Q file")
            continue
        
        user_profile = user_data[0]
        user_components[user_id] = user_profile
        
        # Calculate probabilities
        probs = calculate_ancestry_probabilities(user_profile, country_profiles)
        
        # Bootstrap CI if requested
        ci = {}
        if args.bootstrap > 0:
            ci = bootstrap_confidence_intervals(
                user_profile,
                ref_profiles,
                ref_countries,
                args.bootstrap
            )
        
        # Store results
        for country, prob in sorted(probs.items(), key=lambda x: -x[1]):
            result = {
                'User': user_id,
                'Country': country,
                'ADMIXTURE_Probability': prob * 100
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
    
    # Save user component breakdown
    components_file = args.output.replace('.csv', '_components.csv')
    comp_data = []
    for user_id, profile in user_components.items():
        for i, val in enumerate(profile):
            comp_data.append({
                'User': user_id,
                'Component': f'K{i+1}',
                'Proportion': val,
                'Interpretation': ', '.join(component_interp.get(i+1, ['Unknown']))
            })
    
    pd.DataFrame(comp_data).to_csv(components_file, index=False)
    logger.info(f"Component breakdown saved to: {components_file}")
    
    # Print summary
    print("\n=== ADMIXTURE-Based Ancestry Probabilities ===\n")
    
    for user_id in user_samples:
        user_results = results_df[results_df['User'] == user_id]
        
        print(f"User: {user_id}")
        print("-" * 50)
        
        # Print component breakdown
        print("Ancestral Components:")
        if user_id in user_components:
            for i, val in enumerate(user_components[user_id]):
                interp = ', '.join(component_interp.get(i+1, ['Unknown']))
                print(f"  K{i+1}: {val*100:.1f}% ({interp})")
        
        print("\nCountry Probabilities:")
        for _, row in user_results.head(10).iterrows():
            prob = row['ADMIXTURE_Probability']
            country = row['Country']
            if 'CI_Lower' in row:
                print(f"  {country}: {prob:.1f}% ({row['CI_Lower']:.1f}% - {row['CI_Upper']:.1f}%)")
            else:
                print(f"  {country}: {prob:.1f}%")
        print()


if __name__ == '__main__':
    main()
