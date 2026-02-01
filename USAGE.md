# Usage Guide

This guide provides step-by-step instructions for common usage scenarios.

## Table of Contents

1. [Quick Start](#quick-start)
2. [Understanding Your Results](#understanding-your-results)
3. [Advanced Usage](#advanced-usage)
4. [Troubleshooting](#troubleshooting)
5. [Tips for Better Results](#tips-for-better-results)

## Quick Start

### Step 1: Get Your Raw DNA Data

Download your raw DNA data from:
- **AncestryDNA**: Account Settings → DNA Data → Download DNA Data
- **23andMe**: Account → Raw Data → Download

You'll receive a file named something like `AncestryDNA.txt` or `genome_full.txt`.

### Step 2: Install Dependencies

```bash
pip install -r requirements.txt
```

### Step 3: Run Analysis

```bash
python analyze_ancestry.py path/to/your/AncestryDNA.txt
```

### Step 4: Review Results

The analysis will output a detailed report showing:
- Your estimated European ancestry percentages by country
- Regional summaries
- Uncertainty indicators
- Confidence metrics

## Understanding Your Results

### Ancestry Predictions

Each country is listed with a percentage showing your estimated genetic similarity:

```
British             :  28.3%
Irish               :  22.1%
  ⚠ Moderate uncertainty (±6.2%)
```

- **Percentage**: Your estimated genetic similarity to this population
- **⚠ Warning indicators**: Show when results have high uncertainty
- **Method comparison**: When PCA and ADMIXTURE disagree significantly

### Uncertainty Levels

- **High uncertainty (±10%+)**: Results vary significantly between methods. Common in:
  - Balkans (Serbian, Croatian, Bulgarian overlap heavily)
  - Central/Eastern Europe boundaries
  - Closely related populations

- **Moderate uncertainty (±5-10%)**: Some variation between methods, but generally reliable

- **No warning**: Both methods agree; higher confidence

### Regional Summary

Groups results by major European regions:
- **Western Europe**: British Isles, France, Benelux
- **Southern Europe**: Mediterranean countries
- **Northern Europe**: Scandinavia, Baltics
- **Eastern Europe**: Slavic nations
- **Central Europe**: Germanic nations

### Confidence Metrics

- **0.0-0.5**: Low confidence - results highly uncertain
- **0.5-0.7**: Moderate confidence - interpret with caution
- **0.7-1.0**: Good confidence - methods generally agree

## Advanced Usage

### Adjusting Analysis Parameters

#### More Principal Components (Better Resolution)

```bash
python analyze_ancestry.py your_data.txt --pca-components 15
```

Higher values (10-20) can capture more subtle genetic variation but may introduce noise.

#### More ADMIXTURE Clusters (Finer Granularity)

```bash
python analyze_ancestry.py your_data.txt --admix-clusters 12
```

Higher values (8-15) attempt to detect more ancestral components.

#### Combined Example

```bash
python analyze_ancestry.py your_data.txt \
    --pca-components 15 \
    --admix-clusters 12 \
    --output detailed_report.txt
```

### Using Custom Reference Data

If you have access to European reference datasets:

```bash
python analyze_ancestry.py your_data.txt \
    --reference european_ref.csv \
    --no-synthetic
```

Reference data format (CSV):
```
population,rs1234,rs5678,rs9012,...
British,0,1,2,...
French,1,1,0,...
```

## Troubleshooting

### "Very few SNPs aligned"

**Cause**: Limited overlap between your data and reference SNPs.

**Solutions**:
- Use a larger reference dataset
- Check that your file format is correct
- Ensure your data isn't heavily filtered

### File Format Errors

**Supported formats**:
```
# Tab-delimited
rsid    chromosome    position    genotype
rs12345    1    12345    AG

# Comma-delimited
rsid,chromosome,position,genotype
rs12345,1,12345,AG
```

**Common issues**:
- Missing headers (algorithm will auto-detect)
- Non-standard delimiters (use tab or comma)
- Compressed files (decompress first)

### Low Confidence Results

If confidence < 0.5:
- Results may be highly uncertain
- Try adjusting parameters
- Consider your data may be from non-European ancestry
- Use more reference samples

### High Uncertainty Warnings

**Normal for**:
- Mixed ancestry
- Regions with genetic overlap (Balkans, Central Europe)
- Border populations

**Not normal for**:
- Clearly distinct populations (e.g., Scandinavian vs Mediterranean)
- May indicate data quality issues

## Tips for Better Results

### 1. Understand the Limitations

- This is a DIY tool, not professional testing
- Results are estimates, not definitive proof
- Genetic ancestry ≠ cultural/genealogical ancestry
- European-focused only

### 2. Use Larger Reference Datasets

The synthetic reference data is for demonstration. For better results:
- Use 1000 Genomes European samples
- Academic reference panels
- Minimum 500+ reference samples recommended

### 3. Interpret Regional Results First

If you see:
- 60% Western Europe
- 30% Southern Europe
- 10% Eastern Europe

Focus on these regional patterns before diving into country-level details.

### 4. Cross-Reference with Known Ancestry

If you know you have Spanish ancestry but results show Portuguese:
- These populations are genetically very similar
- Consider regional summaries instead
- Uncertainty warnings help identify these cases

### 5. Multiple Analyses

Run the analysis multiple times with different parameters:
```bash
# Conservative (fewer components)
python analyze_ancestry.py data.txt --pca-components 5 --admix-clusters 5

# Standard
python analyze_ancestry.py data.txt

# Detailed (more components)
python analyze_ancestry.py data.txt --pca-components 15 --admix-clusters 12
```

Compare results - consistent findings are more reliable.

### 6. Understand Genetic Distance

Closely related populations show high overlap:
- British ↔ Irish
- Spanish ↔ Portuguese
- Polish ↔ Czech
- Serbian ↔ Croatian ↔ Bulgarian

This is expected and reflects actual genetic similarity.

## Using the Python API

For integration into your own scripts:

```python
from ancestry_pipeline import (
    GenotypeParser,
    ReferenceDataHandler,
    PCAAnalyzer,
    AdmixtureAnalyzer,
    CountryPredictor
)

# Parse data
parser = GenotypeParser()
user_data = parser.parse_file('genotype.txt')
user_matrix, snps = parser.to_numeric_matrix()

# Generate reference (or load your own)
ref_handler = ReferenceDataHandler()
ref_data, ref_labels = ref_handler.generate_synthetic_reference()

# Align SNPs
user_indices, ref_indices = ref_handler.align_snps(snps)
user_aligned = user_matrix[:, user_indices]
ref_aligned = ref_data[:, ref_indices]

# PCA Analysis
pca = PCAAnalyzer(n_components=10)
pca.fit_reference(ref_aligned, ref_labels)
projection = pca.transform_user(user_aligned)
pca_ancestry = pca.calculate_population_proximity(projection)

# ADMIXTURE Analysis  
admix = AdmixtureAnalyzer(n_clusters=8)
admix.fit_reference(ref_aligned, ref_labels)
admix_ancestry = admix.estimate_ancestry(user_aligned)

# Combine Results
predictor = CountryPredictor()
final = predictor.combine_predictions(pca_ancestry, admix_ancestry)
report = predictor.generate_report(final, ...)
```

## Example Workflows

### Workflow 1: Basic Analysis
```bash
# Download your raw data
# Save as my_dna.txt

# Run analysis
python analyze_ancestry.py my_dna.txt

# Review output
# Look for high-percentage populations
# Check uncertainty warnings
```

### Workflow 2: Detailed Analysis
```bash
# Run with custom parameters
python analyze_ancestry.py my_dna.txt \
    --pca-components 15 \
    --admix-clusters 10 \
    --output detailed_report.txt

# Review saved report
cat detailed_report.txt

# Compare with standard analysis
python analyze_ancestry.py my_dna.txt \
    --output standard_report.txt
    
diff detailed_report.txt standard_report.txt
```

### Workflow 3: Programming Analysis
```bash
# Run example script to understand API
python examples/example_usage.py

# Create your own analysis script
# Use the API for custom analysis
# Integrate with visualization tools
```

## Getting Help

If you encounter issues:

1. Check this guide's [Troubleshooting](#troubleshooting) section
2. Review the main [README.md](README.md)
3. Run the example script: `python examples/example_usage.py`
4. Check that dependencies are installed: `pip install -r requirements.txt`
5. Open an issue on GitHub with:
   - Error message (if any)
   - Command you ran
   - File format sample (first 10 lines)
   - Python version

## Further Reading

- **Principal Component Analysis (PCA)**: Statistical technique for dimensionality reduction
- **ADMIXTURE**: Software for ancestry estimation using maximum likelihood
- **1000 Genomes Project**: Public dataset with European reference populations
- **Population genetics**: Study of genetic variation within and between populations

---

Remember: This tool is for educational purposes. Results should be interpreted as estimates, not definitive proof of ancestry.
