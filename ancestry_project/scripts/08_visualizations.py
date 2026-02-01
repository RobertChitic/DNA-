#!/usr/bin/env python3
"""
Script: 08_visualizations.py
Description: Generate visualizations for ancestry analysis results

Outputs:
- PCA scatter plots (PC1-PC2, PC3-PC4) with user highlighted
- ADMIXTURE stacked barplot
- Side-by-side comparison of all methods
- Interactive HTML components

Usage:
    python 08_visualizations.py --results-dir <dir>
"""

import argparse
import os
import sys
import logging
from typing import Dict, List, Optional

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)

# Color palette for populations
POPULATION_COLORS = {
    'CEU': '#1f77b4',  # Blue
    'GBR': '#ff7f0e',  # Orange
    'FIN': '#2ca02c',  # Green
    'IBS': '#d62728',  # Red
    'TSI': '#9467bd',  # Purple
    'French': '#8c564b',  # Brown
    'Italian': '#e377c2',  # Pink
    'Russian': '#7f7f7f',  # Gray
    'Basque': '#bcbd22',  # Yellow-green
    'Sardinian': '#17becf',  # Cyan
    'Orcadian': '#aec7e8',  # Light blue
    'USER': '#000000',  # Black for user
}

# Country colors
COUNTRY_COLORS = {
    'United Kingdom': '#ff7f0e',
    'United Kingdom (Orkney)': '#aec7e8',
    'Finland': '#2ca02c',
    'Spain': '#d62728',
    'Spain (Basque)': '#bcbd22',
    'Italy': '#9467bd',
    'Italy (Tuscan)': '#e377c2',
    'Italy (Sardinia)': '#17becf',
    'France': '#8c564b',
    'Russia': '#7f7f7f',
    'United States (Utah/Northern European)': '#1f77b4',
}


def load_pca_data(results_dir: str) -> Optional[pd.DataFrame]:
    """Load PCA eigenvec data."""
    eigenvec_file = os.path.join(results_dir, 'pca', 'eigenvec.txt')
    
    if not os.path.exists(eigenvec_file):
        logger.warning(f"Eigenvec file not found: {eigenvec_file}")
        return None
    
    # Try to detect header
    with open(eigenvec_file, 'r') as f:
        first_line = f.readline().strip()
        has_header = 'PC1' in first_line.upper() or 'FID' in first_line.upper()
    
    if has_header:
        df = pd.read_csv(eigenvec_file, sep=r'\s+', header=0)
    else:
        cols = ['FID', 'IID'] + [f'PC{i}' for i in range(1, 25)]
        df = pd.read_csv(eigenvec_file, sep=r'\s+', header=None, names=cols)
    
    df.columns = df.columns.str.upper()
    
    return df


def load_metadata(results_dir: str) -> Optional[pd.DataFrame]:
    """Load sample metadata."""
    metadata_file = os.path.join(os.path.dirname(results_dir), 'reference', 'sample_metadata.tsv')
    
    if not os.path.exists(metadata_file):
        logger.warning(f"Metadata file not found: {metadata_file}")
        return None
    
    df = pd.read_csv(metadata_file, sep='\t')
    df.columns = df.columns.str.lower()
    
    return df


def plot_pca_scatter(
    eigenvec_df: pd.DataFrame,
    metadata_df: Optional[pd.DataFrame],
    output_dir: str,
    user_ids: List[str] = None
):
    """
    Generate PCA scatter plots (PC1-PC2 and PC3-PC4).
    """
    logger.info("Generating PCA scatter plots...")
    
    # Merge with metadata
    if metadata_df is not None:
        merged = eigenvec_df.merge(
            metadata_df,
            left_on='IID',
            right_on='sample_id',
            how='left'
        )
        merged['population'] = merged['population'].fillna('USER')
    else:
        merged = eigenvec_df.copy()
        merged['population'] = 'Unknown'
    
    # Identify user samples
    if user_ids:
        merged.loc[merged['IID'].isin(user_ids), 'population'] = 'USER'
    
    # PC1 vs PC2 plot
    fig, axes = plt.subplots(1, 2, figsize=(14, 6))
    
    # Left plot: PC1 vs PC2
    ax = axes[0]
    for pop in merged['population'].unique():
        if pop == 'USER':
            continue
        mask = merged['population'] == pop
        color = POPULATION_COLORS.get(pop, '#808080')
        ax.scatter(merged[mask]['PC1'], merged[mask]['PC2'], 
                  c=color, label=pop, alpha=0.6, s=30)
    
    # Plot user last (on top)
    user_mask = merged['population'] == 'USER'
    if user_mask.any():
        ax.scatter(merged[user_mask]['PC1'], merged[user_mask]['PC2'],
                  c='black', marker='*', s=300, label='YOU', zorder=5,
                  edgecolors='red', linewidths=2)
    
    ax.set_xlabel('PC1', fontsize=12)
    ax.set_ylabel('PC2', fontsize=12)
    ax.set_title('PCA: PC1 vs PC2', fontsize=14)
    ax.legend(bbox_to_anchor=(1.02, 1), loc='upper left', fontsize=8)
    ax.grid(True, alpha=0.3)
    
    # Right plot: PC3 vs PC4
    ax = axes[1]
    if 'PC3' in merged.columns and 'PC4' in merged.columns:
        for pop in merged['population'].unique():
            if pop == 'USER':
                continue
            mask = merged['population'] == pop
            color = POPULATION_COLORS.get(pop, '#808080')
            ax.scatter(merged[mask]['PC3'], merged[mask]['PC4'],
                      c=color, label=pop, alpha=0.6, s=30)
        
        if user_mask.any():
            ax.scatter(merged[user_mask]['PC3'], merged[user_mask]['PC4'],
                      c='black', marker='*', s=300, label='YOU', zorder=5,
                      edgecolors='red', linewidths=2)
        
        ax.set_xlabel('PC3', fontsize=12)
        ax.set_ylabel('PC4', fontsize=12)
        ax.set_title('PCA: PC3 vs PC4', fontsize=14)
        ax.legend(bbox_to_anchor=(1.02, 1), loc='upper left', fontsize=8)
        ax.grid(True, alpha=0.3)
    
    plt.tight_layout()
    output_file = os.path.join(output_dir, 'pca_plot.png')
    plt.savefig(output_file, dpi=150, bbox_inches='tight')
    plt.close()
    
    logger.info(f"PCA plot saved: {output_file}")


def plot_admixture_barplot(
    q_file: str,
    fam_file: str,
    metadata_df: Optional[pd.DataFrame],
    output_dir: str,
    k: int = None
):
    """
    Generate ADMIXTURE stacked barplot.
    """
    logger.info("Generating ADMIXTURE barplot...")
    
    # Load Q file
    q_data = pd.read_csv(q_file, sep=r'\s+', header=None)
    n_components = q_data.shape[1]
    q_data.columns = [f'K{i+1}' for i in range(n_components)]
    
    # Load FAM file
    fam_data = pd.read_csv(fam_file, sep=r'\s+', header=None,
                          names=['FID', 'IID', 'PID', 'MID', 'Sex', 'Phenotype'])
    
    # Align lengths
    min_len = min(len(fam_data), len(q_data))
    fam_data = fam_data.iloc[:min_len]
    q_data = q_data.iloc[:min_len]
    
    combined = pd.concat([fam_data[['IID']], q_data], axis=1)
    
    # Merge with metadata
    if metadata_df is not None:
        combined = combined.merge(
            metadata_df,
            left_on='IID',
            right_on='sample_id',
            how='left'
        )
        combined['population'] = combined['population'].fillna('USER')
    else:
        combined['population'] = 'Unknown'
    
    # Sort by population, then by dominant component
    combined['dominant'] = q_data.idxmax(axis=1)
    combined = combined.sort_values(['population', 'dominant'])
    
    # Create plot
    k_cols = [col for col in combined.columns if col.startswith('K') and col[1:].isdigit()]
    n_samples = len(combined)
    
    fig, ax = plt.subplots(figsize=(max(12, n_samples * 0.1), 6))
    
    # Color palette for components
    colors = plt.cm.tab20(np.linspace(0, 1, len(k_cols)))
    
    # Stacked bar plot
    bottom = np.zeros(n_samples)
    x = np.arange(n_samples)
    
    for i, k_col in enumerate(k_cols):
        values = combined[k_col].values
        ax.bar(x, values, bottom=bottom, color=colors[i], label=k_col, width=1.0)
        bottom += values
    
    # Add population labels
    pops = combined['population'].values
    pop_changes = [0] + [i+1 for i in range(len(pops)-1) if pops[i] != pops[i+1]] + [len(pops)]
    
    for i in range(len(pop_changes)-1):
        start = pop_changes[i]
        end = pop_changes[i+1]
        mid = (start + end) / 2
        pop_name = pops[start]
        
        # Add vertical line between populations
        if start > 0:
            ax.axvline(x=start-0.5, color='black', linewidth=0.5)
    
    ax.set_xlim(-0.5, n_samples - 0.5)
    ax.set_ylim(0, 1)
    ax.set_xlabel('Individuals (sorted by population)', fontsize=12)
    ax.set_ylabel('Ancestry Proportion', fontsize=12)
    ax.set_title(f'ADMIXTURE Results (K={len(k_cols)})', fontsize=14)
    ax.legend(bbox_to_anchor=(1.02, 1), loc='upper left')
    
    # Remove x-axis ticks (too many individuals)
    ax.set_xticks([])
    
    plt.tight_layout()
    output_file = os.path.join(output_dir, 'admixture_barplot.png')
    plt.savefig(output_file, dpi=150, bbox_inches='tight')
    plt.close()
    
    logger.info(f"ADMIXTURE barplot saved: {output_file}")


def plot_comparison(
    pca_probs_file: str,
    admixture_probs_file: str,
    ensemble_file: str,
    output_dir: str
):
    """
    Generate side-by-side comparison of all methods.
    """
    logger.info("Generating comparison plot...")
    
    # Load probability files
    dfs = {}
    
    if os.path.exists(pca_probs_file):
        pca_df = pd.read_csv(pca_probs_file)
        # Handle different column name conventions
        prob_col = 'PCA_Probability' if 'PCA_Probability' in pca_df.columns else 'PCA_%'
        if prob_col in pca_df.columns:
            dfs['PCA'] = pca_df[['Country', prob_col]].rename(
                columns={prob_col: 'Probability'})
    
    if os.path.exists(admixture_probs_file):
        adm_df = pd.read_csv(admixture_probs_file)
        # Handle different column name conventions
        prob_col = 'ADMIXTURE_Probability' if 'ADMIXTURE_Probability' in adm_df.columns else 'ADMIXTURE_%'
        if prob_col in adm_df.columns:
            dfs['ADMIXTURE'] = adm_df[['Country', prob_col]].rename(
                columns={prob_col: 'Probability'})
    
    if os.path.exists(ensemble_file):
        ens_df = pd.read_csv(ensemble_file)
        # Handle different column name conventions
        prob_col = 'Ensemble_Probability' if 'Ensemble_Probability' in ens_df.columns else 'Ensemble_%'
        if prob_col in ens_df.columns:
            dfs['Ensemble'] = ens_df[['Country', prob_col]].rename(
                columns={prob_col: 'Probability'})
    
    if not dfs:
        logger.warning("No probability files found for comparison plot")
        return
    
    # Combine data
    combined = None
    for method, df in dfs.items():
        df = df.copy()
        df['Method'] = method
        if combined is None:
            combined = df
        else:
            combined = pd.concat([combined, df])
    
    # Get top countries
    top_countries = combined.groupby('Country')['Probability'].max().nlargest(10).index.tolist()
    combined = combined[combined['Country'].isin(top_countries)]
    
    # Create plot
    fig, ax = plt.subplots(figsize=(12, 6))
    
    x = np.arange(len(top_countries))
    width = 0.25
    multiplier = 0
    
    for method, df in dfs.items():
        df_filtered = df[df['Country'].isin(top_countries)]
        df_filtered = df_filtered.set_index('Country').reindex(top_countries).fillna(0)
        
        offset = width * multiplier
        bars = ax.bar(x + offset, df_filtered['Probability'].values, width, label=method)
        multiplier += 1
    
    ax.set_xlabel('Country', fontsize=12)
    ax.set_ylabel('Probability (%)', fontsize=12)
    ax.set_title('Ancestry Estimation Comparison', fontsize=14)
    ax.set_xticks(x + width)
    ax.set_xticklabels(top_countries, rotation=45, ha='right')
    ax.legend()
    ax.grid(True, alpha=0.3, axis='y')
    
    plt.tight_layout()
    output_file = os.path.join(output_dir, 'comparison_plot.png')
    plt.savefig(output_file, dpi=150, bbox_inches='tight')
    plt.close()
    
    logger.info(f"Comparison plot saved: {output_file}")


def main():
    parser = argparse.ArgumentParser(
        description='Generate visualizations for ancestry analysis',
        formatter_class=argparse.RawDescriptionHelpFormatter
    )
    
    parser.add_argument(
        '--results-dir', '-r',
        required=True,
        help='Path to results directory'
    )
    
    parser.add_argument(
        '--user-id',
        help='User sample ID (for highlighting)'
    )
    
    args = parser.parse_args()
    
    results_dir = args.results_dir
    
    # Create output directory
    viz_dir = os.path.join(results_dir, 'visualizations')
    os.makedirs(viz_dir, exist_ok=True)
    
    # Load metadata
    project_dir = os.path.dirname(results_dir)
    metadata_df = load_metadata(results_dir)
    
    # User ID
    user_ids = [args.user_id] if args.user_id else None
    
    # Generate PCA plots
    eigenvec_df = load_pca_data(results_dir)
    if eigenvec_df is not None:
        plot_pca_scatter(eigenvec_df, metadata_df, viz_dir, user_ids)
    
    # Generate ADMIXTURE barplot
    admixture_dir = os.path.join(results_dir, 'admixture')
    optimal_k_file = os.path.join(admixture_dir, 'optimal_k.txt')
    
    if os.path.exists(optimal_k_file):
        with open(optimal_k_file) as f:
            k = int(f.read().strip())
        
        q_file = os.path.join(admixture_dir, f'ld_pruned.{k}.Q')
        fam_file = os.path.join(admixture_dir, 'samples.fam')
        
        if os.path.exists(q_file) and os.path.exists(fam_file):
            plot_admixture_barplot(q_file, fam_file, metadata_df, viz_dir, k)
    
    # Generate comparison plot
    pca_probs = os.path.join(results_dir, 'probabilities', 'pca_probabilities.csv')
    adm_probs = os.path.join(results_dir, 'probabilities', 'admixture_probabilities.csv')
    ensemble = os.path.join(results_dir, 'FINAL_ensemble.csv')
    
    if os.path.exists(pca_probs) or os.path.exists(adm_probs) or os.path.exists(ensemble):
        plot_comparison(pca_probs, adm_probs, ensemble, viz_dir)
    
    # Copy plots to main results directory
    for plot_file in ['pca_plot.png', 'admixture_barplot.png', 'comparison_plot.png']:
        src = os.path.join(viz_dir, plot_file)
        dst = os.path.join(results_dir, plot_file)
        if os.path.exists(src):
            import shutil
            shutil.copy(src, dst)
    
    logger.info("Visualization generation complete!")


if __name__ == '__main__':
    main()
