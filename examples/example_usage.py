#!/usr/bin/env python3
"""
Example script demonstrating ancestry analysis pipeline usage.
"""

from ancestry_pipeline import (
    GenotypeParser,
    ReferenceDataHandler,
    PCAAnalyzer,
    AdmixtureAnalyzer,
    CountryPredictor
)
import numpy as np


def example_with_synthetic_data():
    """
    Example using synthetic reference data.
    This is the recommended approach for testing without real reference datasets.
    """
    print("Example: Ancestry Analysis with Synthetic Reference Data")
    print("=" * 60)
    
    # Create synthetic user data
    print("\n1. Creating synthetic user genotype data...")
    np.random.seed(123)
    user_data = np.random.binomial(2, 0.45, (1, 1000))  # 1 individual, 1000 SNPs
    print(f"   Created data: {user_data.shape[0]} individual, {user_data.shape[1]} SNPs")
    
    # Generate synthetic reference data
    print("\n2. Generating synthetic European reference data...")
    ref_handler = ReferenceDataHandler()
    reference_data, reference_labels = ref_handler.generate_synthetic_reference(
        n_samples_per_pop=20,
        n_snps=1000
    )
    print(f"   Generated: {len(reference_data)} samples, {reference_data.shape[1]} SNPs")
    print(f"   Populations: {len(set(reference_labels))}")
    
    # Run PCA analysis
    print("\n3. Running PCA analysis...")
    pca_analyzer = PCAAnalyzer(n_components=10)
    pca_analyzer.fit_reference(reference_data, np.array(reference_labels))
    
    user_projection = pca_analyzer.transform_user(user_data)
    explained_var = pca_analyzer.get_explained_variance()
    print(f"   Top 3 PCs explain {explained_var[:3].sum()*100:.1f}% of variance")
    
    pca_ancestry = pca_analyzer.calculate_population_proximity(
        user_projection, n_neighbors=50
    )
    print(f"   PCA identified {len(pca_ancestry)} potential ancestries")
    
    # Run ADMIXTURE analysis
    print("\n4. Running ADMIXTURE-like analysis...")
    admixture_analyzer = AdmixtureAnalyzer(n_clusters=8)
    admixture_analyzer.fit_reference(reference_data, np.array(reference_labels))
    
    admixture_ancestry = admixture_analyzer.estimate_ancestry(user_data)
    confidence = admixture_analyzer.estimate_admixture_confidence(user_data)
    print(f"   Confidence: {confidence['overall_confidence']:.2f}")
    
    # Combine results
    print("\n5. Combining results and generating report...")
    predictor = CountryPredictor(uncertainty_threshold=3.0)
    
    combined_predictions = predictor.combine_predictions(
        pca_ancestry, admixture_ancestry, pca_weight=0.5
    )
    
    uncertainty = predictor.calculate_uncertainty(pca_ancestry, admixture_ancestry)
    
    report = predictor.generate_report(
        combined_predictions, uncertainty, confidence
    )
    
    print("\n" + report)
    
    # Show top matches
    print("\nTop 5 Ancestry Matches:")
    top_matches = predictor.get_top_matches(combined_predictions, n=5)
    for i, (country, pct) in enumerate(top_matches, 1):
        print(f"  {i}. {country}: {pct:.1f}%")


def example_parsing_genotype_file():
    """
    Example of parsing a real genotype file.
    """
    print("\nExample: Parsing Genotype File")
    print("=" * 60)
    
    parser = GenotypeParser()
    
    # Try to parse example file if it exists
    try:
        data = parser.parse_file('examples/sample_genotype.txt')
        print(f"✓ Parsed {len(data)} SNPs")
        print(f"\nFirst few SNPs:")
        print(data.head())
        
        print(f"\nChromosome distribution:")
        print(data['chromosome'].value_counts().sort_index())
        
        # Convert to numeric
        numeric_data, snp_ids = parser.to_numeric_matrix()
        print(f"\nNumeric matrix shape: {numeric_data.shape}")
        
    except FileNotFoundError:
        print("Sample genotype file not found.")
        print("Create a file at examples/sample_genotype.txt to test parsing.")


def example_reference_data_stats():
    """
    Example showing reference data statistics.
    """
    print("\nExample: Reference Data Statistics")
    print("=" * 60)
    
    ref_handler = ReferenceDataHandler()
    ref_data, ref_labels = ref_handler.generate_synthetic_reference(
        n_samples_per_pop=20,
        n_snps=500
    )
    
    stats = ref_handler.get_reference_stats()
    print(f"\nReference Dataset Statistics:")
    print(f"  Total samples: {stats['n_samples']}")
    print(f"  Total SNPs: {stats['n_snps']}")
    print(f"  Number of populations: {stats['n_populations']}")
    print(f"\nSamples per population:")
    for pop, count in sorted(stats['populations'].items()):
        print(f"  {pop:20s}: {count:3d} samples")
    
    # Show population centroids
    centroids = ref_handler.get_population_centroids()
    print(f"\nCalculated centroids for {len(centroids)} populations")


if __name__ == '__main__':
    print("ANCESTRY ANALYSIS PIPELINE - USAGE EXAMPLES")
    print("=" * 60)
    
    # Run examples
    example_with_synthetic_data()
    print("\n\n")
    
    example_parsing_genotype_file()
    print("\n\n")
    
    example_reference_data_stats()
    
    print("\n" + "=" * 60)
    print("Examples complete!")
    print("\nFor real analysis, use: python analyze_ancestry.py <genotype_file>")
