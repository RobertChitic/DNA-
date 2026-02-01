#!/usr/bin/env python3
"""
Script: 09_ensemble.py
Description: Combine PCA and ADMIXTURE probabilities with weighted ensemble

Default weights: PCA 30%, ADMIXTURE 70%

Usage:
    python 09_ensemble.py --pca-probs <file> --admixture-probs <file> --output <file>
"""

import argparse
import os
import sys
import logging
from typing import Dict

import pandas as pd
import numpy as np

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)

# Default ensemble weights
DEFAULT_WEIGHTS = {
    'PCA': 0.3,
    'ADMIXTURE': 0.7
}


def load_probabilities(filepath: str, method_name: str) -> pd.DataFrame:
    """
    Load probability file and standardize format.
    """
    logger.info(f"Loading {method_name} probabilities: {filepath}")
    
    if not os.path.exists(filepath):
        logger.warning(f"File not found: {filepath}")
        return pd.DataFrame()
    
    df = pd.read_csv(filepath)
    
    # Standardize column names
    df.columns = df.columns.str.strip()
    
    # Find probability column
    prob_col = None
    for col in df.columns:
        if 'probability' in col.lower():
            prob_col = col
            break
    
    if prob_col is None:
        logger.warning(f"No probability column found in {filepath}")
        return pd.DataFrame()
    
    # Standardize to 'Probability' column
    df = df.rename(columns={prob_col: 'Probability'})
    df['Method'] = method_name
    
    # Select relevant columns
    result = df[['User', 'Country', 'Probability', 'Method']].copy()
    
    logger.info(f"Loaded {len(result)} rows from {method_name}")
    
    return result


def calculate_ensemble(
    pca_df: pd.DataFrame,
    admixture_df: pd.DataFrame,
    weights: Dict[str, float] = None
) -> pd.DataFrame:
    """
    Calculate weighted ensemble of probabilities.
    """
    if weights is None:
        weights = DEFAULT_WEIGHTS
    
    logger.info(f"Ensemble weights: {weights}")
    
    # Normalize weights
    total_weight = sum(weights.values())
    weights = {k: v / total_weight for k, v in weights.items()}
    
    # Get all users and countries
    all_users = set()
    all_countries = set()
    
    for df in [pca_df, admixture_df]:
        if not df.empty:
            all_users.update(df['User'].unique())
            all_countries.update(df['Country'].unique())
    
    # Calculate ensemble for each user
    results = []
    
    for user in all_users:
        user_pca = pca_df[pca_df['User'] == user] if not pca_df.empty else pd.DataFrame()
        user_adm = admixture_df[admixture_df['User'] == user] if not admixture_df.empty else pd.DataFrame()
        
        for country in all_countries:
            pca_prob = 0
            adm_prob = 0
            
            if not user_pca.empty:
                pca_match = user_pca[user_pca['Country'] == country]
                if not pca_match.empty:
                    pca_prob = pca_match['Probability'].values[0]
            
            if not user_adm.empty:
                adm_match = user_adm[user_adm['Country'] == country]
                if not adm_match.empty:
                    adm_prob = adm_match['Probability'].values[0]
            
            # Calculate weighted ensemble
            # Handle cases where one method is missing
            active_weights = {}
            if pca_prob > 0 or not pca_df.empty:
                active_weights['PCA'] = weights['PCA']
            if adm_prob > 0 or not admixture_df.empty:
                active_weights['ADMIXTURE'] = weights['ADMIXTURE']
            
            if not active_weights:
                continue
            
            # Renormalize active weights
            total_active = sum(active_weights.values())
            active_weights = {k: v / total_active for k, v in active_weights.items()}
            
            ensemble_prob = (
                active_weights.get('PCA', 0) * pca_prob +
                active_weights.get('ADMIXTURE', 0) * adm_prob
            )
            
            results.append({
                'User': user,
                'Country': country,
                'PCA_%': pca_prob,
                'ADMIXTURE_%': adm_prob,
                'Ensemble_%': ensemble_prob
            })
    
    results_df = pd.DataFrame(results)
    
    # Sort by ensemble probability
    results_df = results_df.sort_values(['User', 'Ensemble_%'], ascending=[True, False])
    
    return results_df


def main():
    parser = argparse.ArgumentParser(
        description='Combine ancestry probabilities with weighted ensemble',
        formatter_class=argparse.RawDescriptionHelpFormatter
    )
    
    parser.add_argument(
        '--pca-probs', '-p',
        required=True,
        help='Path to PCA probabilities CSV'
    )
    
    parser.add_argument(
        '--admixture-probs', '-a',
        required=True,
        help='Path to ADMIXTURE probabilities CSV'
    )
    
    parser.add_argument(
        '--output', '-o',
        required=True,
        help='Output file for ensemble results'
    )
    
    parser.add_argument(
        '--pca-weight',
        type=float,
        default=DEFAULT_WEIGHTS['PCA'],
        help=f'Weight for PCA method (default: {DEFAULT_WEIGHTS["PCA"]})'
    )
    
    parser.add_argument(
        '--admixture-weight',
        type=float,
        default=DEFAULT_WEIGHTS['ADMIXTURE'],
        help=f'Weight for ADMIXTURE method (default: {DEFAULT_WEIGHTS["ADMIXTURE"]})'
    )
    
    args = parser.parse_args()
    
    # Custom weights
    weights = {
        'PCA': args.pca_weight,
        'ADMIXTURE': args.admixture_weight
    }
    
    # Load probabilities
    pca_df = load_probabilities(args.pca_probs, 'PCA')
    admixture_df = load_probabilities(args.admixture_probs, 'ADMIXTURE')
    
    if pca_df.empty and admixture_df.empty:
        logger.error("No probability data available!")
        sys.exit(1)
    
    # Calculate ensemble
    ensemble_df = calculate_ensemble(pca_df, admixture_df, weights)
    
    # Rename columns for final output
    final_df = ensemble_df.rename(columns={
        'PCA_%': 'PCA_%',
        'ADMIXTURE_%': 'ADMIXTURE_%',
        'Ensemble_%': 'Ensemble_%'
    })
    
    # Also create version with renamed columns for compatibility
    output_df = ensemble_df.rename(columns={
        'Ensemble_%': 'Ensemble_Probability'
    })
    
    # Save results
    output_dir = os.path.dirname(args.output)
    if output_dir:
        os.makedirs(output_dir, exist_ok=True)
    
    final_df.to_csv(args.output, index=False)
    logger.info(f"Results saved to: {args.output}")
    
    # Print summary
    print("\n=== Ensemble Ancestry Results ===\n")
    print(f"Weights: PCA={weights['PCA']*100:.0f}%, ADMIXTURE={weights['ADMIXTURE']*100:.0f}%")
    print()
    
    for user in final_df['User'].unique():
        user_results = final_df[final_df['User'] == user]
        
        print(f"User: {user}")
        print("-" * 60)
        print(f"{'Country':<35} {'PCA':<10} {'ADMIX':<10} {'Ensemble':<10}")
        print("-" * 60)
        
        for _, row in user_results.head(10).iterrows():
            country = row['Country'][:33]
            pca = f"{row['PCA_%']:.1f}%"
            adm = f"{row['ADMIXTURE_%']:.1f}%"
            ens = f"{row['Ensemble_%']:.1f}%"
            print(f"{country:<35} {pca:<10} {adm:<10} {ens:<10}")
        
        print()


if __name__ == '__main__':
    main()
