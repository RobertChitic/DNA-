# European Ancestry Analysis Pipeline

A Python-based DIY ancestry analysis tool that analyzes AncestryDNA raw genotype files against European open reference datasets to produce probabilistic modern country ancestry percentages.

## Overview

This pipeline projects a user's genome into European genetic space using PCA (Principal Component Analysis) and ADMIXTURE-like algorithms, then outputs country-level percentage predictions with uncertainty estimates.

**Important Note**: This is a DIY tool prioritizing exploration and curiosity over perfect accuracy. Results clearly state uncertainty and overlap, especially in regions like the Balkans and Central Europe where genetic boundaries are less distinct.

## Features

- **Genotype File Parsing**: Supports AncestryDNA and 23andMe raw data formats
- **PCA Analysis**: Projects genomes into European genetic space
- **ADMIXTURE Analysis**: Estimates ancestry proportions using clustering algorithms
- **Country-Level Predictions**: Produces probabilistic ancestry percentages for European countries
- **Uncertainty Quantification**: Clearly reports confidence levels and method disagreement
- **Synthetic Reference Data**: Includes synthetic European reference data for testing

## Installation

### Prerequisites

- Python 3.7 or higher
- pip

### Install Dependencies

```bash
pip install -r requirements.txt
```

Required packages:
- numpy
- pandas
- scikit-learn
- matplotlib
- scipy

## Quick Start

### Basic Usage

Analyze a genotype file using synthetic reference data:

```bash
python analyze_ancestry.py examples/sample_genotype.txt
```

### Save Report to File

```bash
python analyze_ancestry.py examples/sample_genotype.txt --output my_ancestry_report.txt
```

### Adjust Analysis Parameters

```bash
python analyze_ancestry.py examples/sample_genotype.txt \
    --pca-components 15 \
    --admix-clusters 10
```

### Using Custom Reference Data

```bash
python analyze_ancestry.py my_genotype.txt \
    --reference path/to/european_reference.csv \
    --no-synthetic
```

## Input File Format

The pipeline accepts raw genotype files from AncestryDNA or 23andMe with the following format:

```
# Comments (optional)
RSID    Chromosome    Position    Genotype
rs4477212    1    82154    AA
rs3094315    1    752566    AG
rs3131972    1    752721    AG
...
```

Supported delimiters: tab or comma

## Output

The pipeline produces a detailed report including:

1. **Ancestry Predictions**: Country-level percentages sorted by likelihood
2. **Uncertainty Indicators**: Warnings for high uncertainty populations
3. **Method Comparison**: PCA vs ADMIXTURE estimates when they differ significantly
4. **Regional Summary**: Grouped results by European regions (Western, Southern, Northern, Eastern, Central)
5. **Confidence Metrics**: Overall confidence score (0-1 scale)

Example output:

```
============================================================
EUROPEAN ANCESTRY ANALYSIS REPORT
============================================================

NOTE: This is a DIY ancestry analysis tool.
Results are estimates and should be interpreted with caution.
Uncertainty and overlap are especially high in regions like
the Balkans and Central Europe.

ANCESTRY PREDICTIONS:
------------------------------------------------------------
British             :  28.3%
Irish               :  22.1%
French              :  18.4%
Spanish             :  12.7%
  ⚠ Moderate uncertainty (±6.2%)
Italian             :   8.9%
German              :   5.8%
Dutch               :   3.8%

REGIONAL SUMMARY:
------------------------------------------------------------
Western Europe              :  78.4%
  - British                 :  28.3%
  - Irish                   :  22.1%
  - French                  :  18.4%
  - German                  :   5.8%
  - Dutch                   :   3.8%
Southern Europe             :  21.6%
  - Spanish                 :  12.7%
  - Italian                 :   8.9%

CONFIDENCE METRICS:
------------------------------------------------------------
Overall Confidence: 0.73 (0-1 scale)
✓ Good confidence in results
```

## Architecture

The pipeline consists of five main modules:

### 1. GenotypeParser (`genotype_parser.py`)
- Parses raw genotype files
- Handles data cleaning and normalization
- Converts genotypes to numeric matrices

### 2. ReferenceDataHandler (`reference_data.py`)
- Manages European reference datasets
- Generates synthetic reference data for testing
- Aligns SNPs between user and reference data

### 3. PCAAnalyzer (`pca_analysis.py`)
- Fits PCA model on reference populations
- Projects user genome into genetic space
- Calculates population proximity scores

### 4. AdmixtureAnalyzer (`admixture_analysis.py`)
- Implements ADMIXTURE-like clustering
- Estimates ancestry proportions
- Provides confidence metrics

### 5. CountryPredictor (`country_predictor.py`)
- Combines PCA and ADMIXTURE results
- Calculates uncertainty metrics
- Generates human-readable reports

## Example Usage (Python API)

```python
from ancestry_pipeline import (
    GenotypeParser,
    ReferenceDataHandler,
    PCAAnalyzer,
    AdmixtureAnalyzer,
    CountryPredictor
)
import numpy as np

# Parse genotype file
parser = GenotypeParser()
user_data = parser.parse_file('my_genotype.txt')
user_numeric, snp_ids = parser.to_numeric_matrix()

# Load reference data
ref_handler = ReferenceDataHandler()
ref_data, ref_labels = ref_handler.generate_synthetic_reference(
    n_samples_per_pop=20,
    n_snps=1000
)

# Run PCA analysis
pca_analyzer = PCAAnalyzer(n_components=10)
pca_analyzer.fit_reference(ref_data, np.array(ref_labels))
user_projection = pca_analyzer.transform_user(user_numeric)
pca_ancestry = pca_analyzer.calculate_population_proximity(user_projection)

# Run ADMIXTURE analysis
admixture_analyzer = AdmixtureAnalyzer(n_clusters=8)
admixture_analyzer.fit_reference(ref_data, np.array(ref_labels))
admixture_ancestry = admixture_analyzer.estimate_ancestry(user_numeric)

# Combine and report
predictor = CountryPredictor()
combined = predictor.combine_predictions(pca_ancestry, admixture_ancestry)
uncertainty = predictor.calculate_uncertainty(pca_ancestry, admixture_ancestry)
report = predictor.generate_report(combined, uncertainty)
print(report)
```

## Limitations and Caveats

1. **Not Medical Grade**: This is a DIY tool for educational/curiosity purposes only
2. **Synthetic Reference Data**: Default reference data is synthetic; real datasets improve accuracy
3. **European Focus**: Optimized for European ancestry only
4. **SNP Coverage**: Results depend on SNP overlap between user and reference data
5. **Regional Ambiguity**: High uncertainty in genetically similar regions (e.g., Balkans, Central Europe)
6. **Simplified Algorithms**: Uses simplified versions of PCA/ADMIXTURE for accessibility

## Using Real Reference Data

For better accuracy, use real European reference datasets:

1. **1000 Genomes Project**: Public dataset with European populations
2. **Human Genome Diversity Project (HGDP)**: Includes European samples
3. **European reference panels**: Various academic sources

Reference data should be in CSV format:
```
population,rs1,rs2,rs3,...
British,0,1,2,...
French,1,1,0,...
...
```

## Testing

Run the example script to verify installation:

```bash
python examples/example_usage.py
```

This will:
- Demonstrate synthetic data generation
- Parse the sample genotype file
- Show reference data statistics
- Run a complete analysis pipeline

## Contributing

Contributions are welcome! Areas for improvement:

- Support for additional reference data formats (VCF, PLINK)
- Integration with public reference datasets
- Visualization of PCA projections
- Web interface
- Mobile app support
- Non-European ancestry support

## License

This project is provided as-is for educational purposes.

## Acknowledgments

- Inspired by professional ancestry testing services
- Uses standard bioinformatics approaches (PCA, ADMIXTURE)
- Built with scientific Python stack (NumPy, pandas, scikit-learn)

## Disclaimer

This tool is for educational and entertainment purposes only. Results should not be used for medical decisions, legal purposes, or as definitive proof of ancestry. For accurate ancestry information, consider professional genetic testing services.

## Contact

For questions, issues, or suggestions, please open an issue on GitHub.

---

**Remember**: Genetics is complex, and ancestry is only one small part of your identity. This tool provides estimates based on current genetic data and should be interpreted with appropriate skepticism and curiosity.