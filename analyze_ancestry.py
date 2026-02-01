#!/usr/bin/env python3
"""
Main pipeline script for European ancestry analysis.
"""

import argparse
import sys
import numpy as np
from pathlib import Path

from ancestry_pipeline import (
    GenotypeParser,
    ReferenceDataHandler,
    PCAAnalyzer,
    AdmixtureAnalyzer,
    CountryPredictor
)


def run_ancestry_analysis(genotype_file: str, 
                          reference_file: str = None,
                          use_synthetic: bool = True,
                          n_components: int = 10,
                          n_clusters: int = 8,
                          output_file: str = None):
    """
    Run complete ancestry analysis pipeline.
    
    Args:
        genotype_file: Path to user's genotype file
        reference_file: Path to reference dataset (optional)
        use_synthetic: Whether to use synthetic reference data
        n_components: Number of PCA components
        n_clusters: Number of ADMIXTURE clusters
        output_file: Path to save report (optional)
    """
    print("=" * 60)
    print("EUROPEAN ANCESTRY ANALYSIS PIPELINE")
    print("=" * 60)
    print()
    
    # Step 1: Parse genotype file
    print("Step 1: Parsing genotype file...")
    parser = GenotypeParser()
    
    try:
        user_data_df = parser.parse_file(genotype_file)
        print(f"✓ Loaded {len(user_data_df)} SNPs from genotype file")
    except Exception as e:
        print(f"✗ Error parsing genotype file: {e}")
        return None
    
    # Step 2: Load or generate reference data
    print("\nStep 2: Loading reference data...")
    ref_handler = ReferenceDataHandler()
    
    if use_synthetic or reference_file is None:
        print("  Using synthetic reference data (for demonstration)...")
        reference_data, reference_labels = ref_handler.generate_synthetic_reference(
            n_samples_per_pop=20,
            n_snps=1000
        )
        print(f"✓ Generated synthetic reference: {len(reference_data)} samples, "
              f"{reference_data.shape[1]} SNPs")
    else:
        try:
            reference_data, reference_labels = ref_handler.load_reference_data(
                reference_file
            )
            print(f"✓ Loaded reference data: {len(reference_data)} samples")
        except Exception as e:
            print(f"✗ Error loading reference data: {e}")
            return None
    
    stats = ref_handler.get_reference_stats()
    print(f"  Populations: {stats['n_populations']}")
    
    # Step 3: Align SNPs
    print("\nStep 3: Aligning SNPs between user and reference...")
    user_numeric, user_snps = parser.to_numeric_matrix()
    
    # For synthetic data, we need to simulate alignment
    # In real scenario, would align on RSIDs
    if use_synthetic:
        # Trim user data to match reference size for demo
        n_snps_to_use = min(len(user_snps), reference_data.shape[1])
        user_numeric = user_numeric[:, :n_snps_to_use]
        reference_data = reference_data[:, :n_snps_to_use]
        print(f"✓ Aligned {n_snps_to_use} SNPs")
    else:
        user_indices, ref_indices = ref_handler.align_snps(user_snps)
        user_numeric = user_numeric[:, user_indices]
        reference_data = reference_data[:, ref_indices]
        print(f"✓ Aligned {len(user_indices)} overlapping SNPs")
    
    if user_numeric.shape[1] < 100:
        print("⚠ Warning: Very few SNPs aligned. Results may be unreliable.")
    
    # Step 4: PCA Analysis
    print("\nStep 4: Running PCA analysis...")
    pca_analyzer = PCAAnalyzer(n_components=n_components)
    pca_analyzer.fit_reference(reference_data, np.array(reference_labels))
    
    user_projection = pca_analyzer.transform_user(user_numeric)
    explained_var = pca_analyzer.get_explained_variance()
    print(f"✓ PCA complete. Top 3 components explain "
          f"{explained_var[:3].sum()*100:.1f}% of variance")
    
    # Get PCA-based ancestry
    pca_ancestry = pca_analyzer.calculate_population_proximity(
        user_projection,
        n_neighbors=50
    )
    print(f"  PCA identified {len(pca_ancestry)} potential ancestries")
    
    # Step 5: ADMIXTURE Analysis
    print("\nStep 5: Running ADMIXTURE-like analysis...")
    admixture_analyzer = AdmixtureAnalyzer(n_clusters=n_clusters)
    admixture_analyzer.fit_reference(reference_data, np.array(reference_labels))
    
    admixture_ancestry = admixture_analyzer.estimate_ancestry(user_numeric)
    confidence = admixture_analyzer.estimate_admixture_confidence(user_numeric)
    print(f"✓ ADMIXTURE complete. Confidence: {confidence['overall_confidence']:.2f}")
    print(f"  Identified {len(admixture_ancestry)} ancestry components")
    
    # Step 6: Combine results and generate report
    print("\nStep 6: Combining results and generating report...")
    predictor = CountryPredictor(uncertainty_threshold=3.0)
    
    combined_predictions = predictor.combine_predictions(
        pca_ancestry,
        admixture_ancestry,
        pca_weight=0.5
    )
    
    uncertainty = predictor.calculate_uncertainty(pca_ancestry, admixture_ancestry)
    
    report = predictor.generate_report(
        combined_predictions,
        uncertainty,
        confidence
    )
    
    print("\n")
    print(report)
    
    # Save report if requested
    if output_file:
        with open(output_file, 'w') as f:
            f.write(report)
        print(f"\n✓ Report saved to {output_file}")
    
    return {
        'combined_predictions': combined_predictions,
        'uncertainty': uncertainty,
        'confidence': confidence,
        'report': report
    }


def main():
    """Main entry point."""
    parser = argparse.ArgumentParser(
        description='Analyze AncestryDNA files against European reference datasets',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Analyze a genotype file with synthetic reference data
  python analyze_ancestry.py my_genotype.txt
  
  # Analyze with custom reference data
  python analyze_ancestry.py my_genotype.txt --reference euro_reference.csv
  
  # Save report to file
  python analyze_ancestry.py my_genotype.txt --output my_report.txt
  
  # Adjust analysis parameters
  python analyze_ancestry.py my_genotype.txt --pca-components 15 --admix-clusters 10
        """
    )
    
    parser.add_argument(
        'genotype_file',
        help='Path to AncestryDNA or 23andMe raw genotype file'
    )
    
    parser.add_argument(
        '--reference',
        help='Path to reference dataset (CSV format)',
        default=None
    )
    
    parser.add_argument(
        '--no-synthetic',
        action='store_true',
        help='Do not use synthetic data (requires --reference)'
    )
    
    parser.add_argument(
        '--pca-components',
        type=int,
        default=10,
        help='Number of PCA components to compute (default: 10)'
    )
    
    parser.add_argument(
        '--admix-clusters',
        type=int,
        default=8,
        help='Number of ADMIXTURE clusters (default: 8)'
    )
    
    parser.add_argument(
        '--output',
        help='Path to save output report'
    )
    
    args = parser.parse_args()
    
    # Validate arguments
    if args.no_synthetic and not args.reference:
        parser.error("--no-synthetic requires --reference to be specified")
    
    if not Path(args.genotype_file).exists():
        print(f"Error: Genotype file not found: {args.genotype_file}")
        sys.exit(1)
    
    # Run analysis
    try:
        results = run_ancestry_analysis(
            genotype_file=args.genotype_file,
            reference_file=args.reference,
            use_synthetic=not args.no_synthetic,
            n_components=args.pca_components,
            n_clusters=args.admix_clusters,
            output_file=args.output
        )
        
        if results is None:
            sys.exit(1)
            
    except Exception as e:
        print(f"\n✗ Analysis failed: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)


if __name__ == '__main__':
    main()
