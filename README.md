# DNA- 🧬

DIY Ancestry Analysis Pipeline for analyzing AncestryDNA raw genotype files against European reference populations.

## Quick Start

See the full documentation in [`ancestry_project/README.md`](ancestry_project/README.md).

```bash
cd ancestry_project
conda env create -f environment.yml
conda activate ancestry_analysis
./run_all.sh /path/to/your/AncestryDNA.txt YourName
```

## Project Structure

```
ancestry_project/
├── data/
│   ├── raw/          # User uploads, reference datasets
│   └── processed/    # PLINK binary files
├── reference/        # Population mappings, metadata
├── results/          # PCA, ADMIXTURE, predictions
├── scripts/          # Analysis scripts (01-10)
├── environment.yml   # Conda dependencies
├── run_all.sh        # Master pipeline script
└── README.md         # Detailed documentation
```

## License

MIT License