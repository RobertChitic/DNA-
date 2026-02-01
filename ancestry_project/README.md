# 🧬 DIY Ancestry Analysis Pipeline

A comprehensive Python-based pipeline for analyzing AncestryDNA raw genotype files against European reference populations. This tool predicts modern country ancestry percentages using PCA and ADMIXTURE methods, providing more detailed results than commercial ancestry services.

## ⚠️ Important Disclaimer

**These results are statistical estimates, not definitive ancestry percentages.**

Genetic ancestry analysis has inherent limitations:
- Reference panels represent only a subset of genetic diversity within each country
- Modern country borders do not reflect historical population movements
- Different reference datasets and methods may produce different results
- These estimates reflect genetic similarity to reference populations, not documented genealogy

Use these results for educational and entertainment purposes only. For medical or legal purposes, consult qualified professionals.

---

## 📋 Table of Contents

- [Features](#features)
- [Requirements](#requirements)
- [Installation](#installation)
- [Quick Start](#quick-start)
- [Pipeline Overview](#pipeline-overview)
- [Configuration](#configuration)
- [Output Files](#output-files)
- [Interpreting Results](#interpreting-results)
- [Advanced Usage](#advanced-usage)
- [Troubleshooting](#troubleshooting)
- [References](#references)

---

## ✨ Features

- **Dual Analysis Methods**: Combines PCA (Principal Component Analysis) and ADMIXTURE for robust ancestry estimation
- **European Reference Populations**: 
  - 1000 Genomes: CEU, GBR, FIN, IBS, TSI
  - HGDP: French, Italian, Russian, Basque, Sardinian, Orcadian
- **Ensemble Predictions**: Weighted combination of methods (default: 30% PCA, 70% ADMIXTURE)
- **Bootstrap Confidence Intervals**: Optional resampling for uncertainty quantification
- **Interactive Visualizations**: HTML dashboard with Plotly charts
- **Automatic Strand Flip Handling**: Robust data processing for varied input formats

---

## 💻 Requirements

### System Requirements

| Resource | Minimum | Recommended |
|----------|---------|-------------|
| RAM | 8 GB | 16+ GB |
| Disk Space | 75 GB | 100+ GB |
| CPU | 4 cores | 8+ cores |
| Runtime | 4-6 hours | 2-4 hours |

### Software Dependencies

- **Python** 3.10+
- **PLINK** 1.9
- **ADMIXTURE** 1.3.0
- **bcftools** (for data preparation)
- **wget** (for data download)

### Python Packages

```
pandas >= 2.0
numpy >= 1.24
scipy >= 1.11
scikit-learn >= 1.3
matplotlib >= 3.7
seaborn >= 0.12
plotly >= 5.18
tqdm >= 4.65
```

---

## 🔧 Installation

### Option 1: Using Conda (Recommended)

```bash
# Clone the repository
git clone https://github.com/yourusername/ancestry_project.git
cd ancestry_project

# Create conda environment
conda env create -f environment.yml

# Activate environment
conda activate ancestry_analysis

# Install ADMIXTURE (if not included)
wget http://dalexander.github.io/admixture/binaries/admixture_linux-1.3.0.tar.gz
tar xzf admixture_linux-1.3.0.tar.gz
sudo mv admixture /usr/local/bin/
```

### Option 2: Manual Installation

```bash
# Install system dependencies (Ubuntu/Debian)
sudo apt-get update
sudo apt-get install plink1.9 bcftools wget

# Create virtual environment
python3 -m venv venv
source venv/bin/activate

# Install Python packages
pip install pandas numpy scipy scikit-learn matplotlib seaborn plotly tqdm kaleido

# Install ADMIXTURE
wget http://dalexander.github.io/admixture/binaries/admixture_linux-1.3.0.tar.gz
tar xzf admixture_linux-1.3.0.tar.gz
sudo mv admixture /usr/local/bin/
```

### Verify Installation

```bash
# Check all tools are available
plink --version
admixture
python3 -c "import pandas; import plotly; print('Python packages OK')"
```

---

## 🚀 Quick Start

### 1. Download Your Raw DNA Data

Export your raw DNA data from AncestryDNA:
1. Log into your Ancestry.com account
2. Go to DNA > Settings > Download Raw DNA Data
3. Save the file (typically `AncestryDNA.txt`)

### 2. Run the Analysis

```bash
# Navigate to project directory
cd ancestry_project

# Make scripts executable
chmod +x run_all.sh scripts/*.sh

# Run complete pipeline
./run_all.sh /path/to/your/AncestryDNA.txt YourName
```

### 3. View Results

Open `results/interactive_dashboard.html` in your web browser.

---

## 📊 Pipeline Overview

The pipeline consists of 10 sequential steps:

```
┌─────────────────────────────────────────────────────────────────┐
│  Step 1: Download Reference Data (01_download_data.sh)          │
│  - 1000 Genomes Phase 3 European populations                    │
│  - HGDP European populations                                    │
└──────────────────────────────┬──────────────────────────────────┘
                               ▼
┌─────────────────────────────────────────────────────────────────┐
│  Step 2: Convert AncestryDNA File (02_convert_ancestry_file.py) │
│  - Parse raw genotype file                                      │
│  - Generate PLINK binary format (.bed/.bim/.fam)                │
└──────────────────────────────┬──────────────────────────────────┘
                               ▼
┌─────────────────────────────────────────────────────────────────┐
│  Step 3: Merge Datasets (03_merge_references.sh)                │
│  - Align user data with references                              │
│  - Handle strand flips automatically                            │
└──────────────────────────────┬──────────────────────────────────┘
                               ▼
┌─────────────────────────────────────────────────────────────────┐
│  Step 4: PCA Analysis (04_run_pca.sh)                           │
│  - LD pruning (indep-pairwise 50 5 0.2)                         │
│  - Compute 20 principal components                              │
└──────────────────────────────┬──────────────────────────────────┘
                               ▼
┌─────────────────────────────────────────────────────────────────┐
│  Step 5: ADMIXTURE Analysis (05_run_admixture.sh)               │
│  - Run K=2 to K=15                                              │
│  - Select optimal K via cross-validation                        │
└──────────────────────────────┬──────────────────────────────────┘
                               ▼
┌─────────────────────────────────────────────────────────────────┐
│  Step 6: PCA Probabilities (06_pca_probabilities.py)            │
│  - Calculate Euclidean distances in PC space                    │
│  - Convert to probabilities via Gaussian kernel                 │
└──────────────────────────────┬──────────────────────────────────┘
                               ▼
┌─────────────────────────────────────────────────────────────────┐
│  Step 7: ADMIXTURE Probabilities (07_admixture_probabilities.py)│
│  - Calculate cosine similarity with reference profiles          │
│  - Convert to country-level probabilities                       │
└──────────────────────────────┬──────────────────────────────────┘
                               ▼
┌─────────────────────────────────────────────────────────────────┐
│  Step 8: Ensemble Combination (09_ensemble.py)                  │
│  - Weighted average: 30% PCA + 70% ADMIXTURE                    │
└──────────────────────────────┬──────────────────────────────────┘
                               ▼
┌─────────────────────────────────────────────────────────────────┐
│  Step 9-10: Visualizations & Report                             │
│  - PCA scatter plots, ADMIXTURE barplots                        │
│  - Interactive HTML dashboard                                   │
└─────────────────────────────────────────────────────────────────┘
```

---

## ⚙️ Configuration

### Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `N_PCS` | 10 | Number of principal components to use |
| `BANDWIDTH` | 0.5 | Gaussian kernel bandwidth for PCA probabilities |
| `MIN_K` | 2 | Minimum K value for ADMIXTURE |
| `MAX_K` | 15 | Maximum K value for ADMIXTURE |
| `BOOTSTRAP` | 0 | Bootstrap iterations (0 = disabled) |
| `PCA_WEIGHT` | 0.3 | Weight for PCA in ensemble |
| `ADMIXTURE_WEIGHT` | 0.7 | Weight for ADMIXTURE in ensemble |

### Setting Parameters

```bash
# Using environment variables
export N_PCS=15
export BANDWIDTH=0.3  # Lower bandwidth = narrower kernel = more localized comparisons
export BOOTSTRAP=1000
./run_all.sh your_data.txt

# Or inline
N_PCS=15 BOOTSTRAP=1000 ./run_all.sh your_data.txt
```

---

## 📁 Output Files

After running the pipeline, you'll find the following outputs:

### Main Results

| File | Description |
|------|-------------|
| `results/FINAL_ensemble.csv` | Final ancestry percentages with all methods |
| `results/interactive_dashboard.html` | Interactive visualization dashboard |

### Visualizations

| File | Description |
|------|-------------|
| `results/pca_plot.png` | PCA scatter plot with your position |
| `results/admixture_barplot.png` | ADMIXTURE component visualization |
| `results/comparison_plot.png` | Method comparison chart |

### Intermediate Results

| Directory | Contents |
|-----------|----------|
| `results/pca/` | PCA eigenvectors, eigenvalues, scree plot |
| `results/admixture/` | Q files, P files, CV error plot |
| `results/probabilities/` | Per-method probability files |

---

## 📖 Interpreting Results

### Understanding the Output

The `FINAL_ensemble.csv` file contains:

```csv
User,Country,PCA_%,ADMIXTURE_%,Ensemble_%
USER001,United Kingdom,25.5,28.2,27.4
USER001,Finland,18.3,15.1,16.1
USER001,Italy,12.4,14.8,14.1
...
```

- **PCA_%**: Probability based on position in principal component space
- **ADMIXTURE_%**: Probability based on ancestral component similarity
- **Ensemble_%**: Weighted combination (recommended for interpretation)

### PCA Plot Interpretation

- Your sample is marked with a star (★)
- Closer proximity to a population cluster suggests higher genetic similarity
- PC1 typically separates North-South European populations
- PC2 often separates East-West populations

### ADMIXTURE Interpretation

- Each component (K1, K2, ...) represents an inferred ancestral population
- Components are not real historical populations but mathematical constructs
- The optimal K is chosen by minimizing cross-validation error

### Confidence Intervals

If `BOOTSTRAP > 0`, results include 95% confidence intervals:
- Wide intervals indicate uncertainty in the estimate
- Overlapping intervals between countries mean those estimates are not significantly different

---

## 🔬 Advanced Usage

### Running Individual Steps

```bash
# Run only PCA analysis
./scripts/04_run_pca.sh data/processed/merged_all 20

# Run only ADMIXTURE
./scripts/05_run_admixture.sh results/pca/ld_pruned 2 12

# Regenerate visualizations
python scripts/08_visualizations.py --results-dir results
```

### Custom Reference Panels

To add custom reference samples:

1. Prepare your reference data in PLINK format
2. Create a metadata file with columns: `sample_id`, `population`, `country`
3. Modify `03_merge_references.sh` to include your reference files

### Adjusting Ensemble Weights

Based on your data quality and sample size:

```bash
# More weight on ADMIXTURE (default)
PCA_WEIGHT=0.3 ADMIXTURE_WEIGHT=0.7 ./run_all.sh data.txt

# Equal weights
PCA_WEIGHT=0.5 ADMIXTURE_WEIGHT=0.5 ./run_all.sh data.txt

# More weight on PCA (if ADMIXTURE has convergence issues)
PCA_WEIGHT=0.6 ADMIXTURE_WEIGHT=0.4 ./run_all.sh data.txt
```

### Re-running Analysis

```bash
# Clean all intermediate files and re-run
./run_all.sh --clean your_data.txt

# Or manually remove step markers
rm results/.step*_complete
./run_all.sh your_data.txt
```

---

## 🔧 Troubleshooting

### Common Issues

#### "PLINK not found"
```bash
# Install via conda
conda install -c bioconda plink

# Or download manually
wget https://s3.amazonaws.com/plink1-assets/plink_linux_x86_64_20231211.zip
unzip plink_linux_x86_64_20231211.zip
sudo mv plink /usr/local/bin/
```

#### "ADMIXTURE not found"
```bash
wget http://dalexander.github.io/admixture/binaries/admixture_linux-1.3.0.tar.gz
tar xzf admixture_linux-1.3.0.tar.gz
sudo mv admixture /usr/local/bin/
```

#### "Out of memory"
- Reduce `MAX_K` to 10
- Run on a machine with more RAM
- Process chromosomes separately

#### "Strand flip errors"
The pipeline handles most strand flips automatically. If issues persist:
- Check input file format matches expected AncestryDNA format
- Ensure alleles are in standard notation (A/C/G/T)

### Error Messages

| Error | Cause | Solution |
|-------|-------|----------|
| "No common SNPs" | Input format mismatch | Check input file header format |
| "CV error not converging" | Too few samples | Reduce MAX_K |
| "Memory allocation failed" | Insufficient RAM | Use machine with more memory |

---

## 📚 References

### Methods

1. **PCA**: Price et al. (2006). "Principal components analysis corrects for stratification in genome-wide association studies." *Nature Genetics*.

2. **ADMIXTURE**: Alexander et al. (2009). "Fast model-based estimation of ancestry in unrelated individuals." *Genome Research*.

### Reference Datasets

1. **1000 Genomes Project**: The 1000 Genomes Project Consortium (2015). "A global reference for human genetic variation." *Nature*.

2. **HGDP**: Bergström et al. (2020). "Insights into human genetic variation and population history from 929 diverse genomes." *Science*.

### Software

- PLINK: https://www.cog-genomics.org/plink/
- ADMIXTURE: http://dalexander.github.io/admixture/
- bcftools: http://samtools.github.io/bcftools/

---

## 📝 License

This project is released under the MIT License. See [LICENSE](LICENSE) for details.

---

## 🤝 Contributing

Contributions are welcome! Please feel free to submit a Pull Request.

---

## 📬 Support

For questions or issues:
1. Check the [Troubleshooting](#troubleshooting) section
2. Open an issue on GitHub
3. Review existing issues for similar problems

---

**Note**: This tool is for educational and research purposes. Results should not be used for medical, legal, or immigration purposes without consulting appropriate professionals.
