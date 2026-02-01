#!/bin/bash
#===============================================================================
# Script: 04_run_pca.sh
# Description: Run PCA analysis using PLINK with LD pruning
#
# Usage: ./04_run_pca.sh [input_prefix] [n_pcs]
#
# Input:
#   - Merged PLINK files from 03_merge_references.sh
#
# Output:
#   - results/pca/merged_all.eigenvec (PC coordinates)
#   - results/pca/merged_all.eigenval (eigenvalues)
#   - results/pca/scree_plot.png (variance explained)
#===============================================================================

set -euo pipefail

# Configuration
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_DIR="$(dirname "$SCRIPT_DIR")"
PROCESSED_DIR="${PROJECT_DIR}/data/processed"
RESULTS_DIR="${PROJECT_DIR}/results/pca"
LOG_FILE="${PROJECT_DIR}/results/04_pca.log"

# Parameters
INPUT_PREFIX="${1:-${PROCESSED_DIR}/merged_all}"
N_PCS="${2:-20}"

# LD pruning parameters
LD_WINDOW=50
LD_STEP=5
LD_R2=0.2

# Create directories
mkdir -p "$RESULTS_DIR" "$(dirname "$LOG_FILE")"

# Initialize log
echo "============================================" | tee "$LOG_FILE"
echo "Ancestry Analysis - PCA" | tee -a "$LOG_FILE"
echo "Started: $(date)" | tee -a "$LOG_FILE"
echo "============================================" | tee -a "$LOG_FILE"
echo "Input: ${INPUT_PREFIX}" | tee -a "$LOG_FILE"
echo "Number of PCs: ${N_PCS}" | tee -a "$LOG_FILE"
echo "LD pruning: window=${LD_WINDOW}, step=${LD_STEP}, r2=${LD_R2}" | tee -a "$LOG_FILE"

#-------------------------------------------------------------------------------
# Check input files
#-------------------------------------------------------------------------------
check_inputs() {
    if [[ ! -f "${INPUT_PREFIX}.bed" ]] || [[ ! -f "${INPUT_PREFIX}.bim" ]] || [[ ! -f "${INPUT_PREFIX}.fam" ]]; then
        echo "ERROR: Input PLINK files not found: ${INPUT_PREFIX}" | tee -a "$LOG_FILE"
        exit 1
    fi
    
    local n_variants=$(wc -l < "${INPUT_PREFIX}.bim")
    local n_samples=$(wc -l < "${INPUT_PREFIX}.fam")
    echo "Input dataset: ${n_variants} variants, ${n_samples} samples" | tee -a "$LOG_FILE"
}

#-------------------------------------------------------------------------------
# Quality control filtering
#-------------------------------------------------------------------------------
qc_filter() {
    echo "" | tee -a "$LOG_FILE"
    echo "=== Quality Control Filtering ===" | tee -a "$LOG_FILE"
    
    local qc_prefix="${RESULTS_DIR}/qc_filtered"
    
    if [[ ! -f "${qc_prefix}.bed" ]]; then
        plink --bfile "$INPUT_PREFIX" \
            --maf 0.01 \
            --geno 0.05 \
            --mind 0.1 \
            --hwe 1e-6 \
            --make-bed \
            --out "$qc_prefix" \
            --allow-no-sex \
            >> "$LOG_FILE" 2>&1
    else
        echo "QC filtered files already exist, skipping." | tee -a "$LOG_FILE"
    fi
    
    local n_variants=$(wc -l < "${qc_prefix}.bim")
    local n_samples=$(wc -l < "${qc_prefix}.fam")
    echo "After QC: ${n_variants} variants, ${n_samples} samples" | tee -a "$LOG_FILE"
    
    echo "$qc_prefix"
}

#-------------------------------------------------------------------------------
# LD pruning
#-------------------------------------------------------------------------------
ld_prune() {
    local input_prefix="$1"
    
    echo "" | tee -a "$LOG_FILE"
    echo "=== LD Pruning ===" | tee -a "$LOG_FILE"
    
    local pruned_prefix="${RESULTS_DIR}/ld_pruned"
    
    if [[ ! -f "${pruned_prefix}.bed" ]]; then
        # Calculate LD and identify SNPs to keep
        echo "Calculating LD structure..." | tee -a "$LOG_FILE"
        plink --bfile "$input_prefix" \
            --indep-pairwise $LD_WINDOW $LD_STEP $LD_R2 \
            --out "${RESULTS_DIR}/ld_prune" \
            --allow-no-sex \
            >> "$LOG_FILE" 2>&1
        
        local n_pruned=$(wc -l < "${RESULTS_DIR}/ld_prune.prune.in")
        echo "SNPs after LD pruning: ${n_pruned}" | tee -a "$LOG_FILE"
        
        # Extract pruned SNPs
        plink --bfile "$input_prefix" \
            --extract "${RESULTS_DIR}/ld_prune.prune.in" \
            --make-bed \
            --out "$pruned_prefix" \
            --allow-no-sex \
            >> "$LOG_FILE" 2>&1
    else
        echo "LD pruned files already exist, skipping." | tee -a "$LOG_FILE"
    fi
    
    echo "$pruned_prefix"
}

#-------------------------------------------------------------------------------
# Run PCA
#-------------------------------------------------------------------------------
run_pca() {
    local input_prefix="$1"
    
    echo "" | tee -a "$LOG_FILE"
    echo "=== Running PCA ===" | tee -a "$LOG_FILE"
    
    local pca_output="${RESULTS_DIR}/pca_results"
    
    if [[ ! -f "${pca_output}.eigenvec" ]]; then
        plink --bfile "$input_prefix" \
            --pca $N_PCS header tabs \
            --out "$pca_output" \
            --allow-no-sex \
            >> "$LOG_FILE" 2>&1
    else
        echo "PCA results already exist, skipping." | tee -a "$LOG_FILE"
    fi
    
    echo "PCA completed" | tee -a "$LOG_FILE"
    
    # Copy results to final location
    cp "${pca_output}.eigenvec" "${RESULTS_DIR}/eigenvec.txt"
    cp "${pca_output}.eigenval" "${RESULTS_DIR}/eigenval.txt"
}

#-------------------------------------------------------------------------------
# Generate scree plot
#-------------------------------------------------------------------------------
generate_scree_plot() {
    echo "" | tee -a "$LOG_FILE"
    echo "=== Generating Scree Plot ===" | tee -a "$LOG_FILE"
    
    python3 << 'EOF'
import matplotlib.pyplot as plt
import numpy as np
import os

# Read eigenvalues
script_dir = os.path.dirname(os.path.abspath('$0'))
results_dir = os.path.join(os.path.dirname(script_dir), 'results', 'pca')
eigenval_file = os.path.join(results_dir, 'eigenval.txt')

eigenvalues = np.loadtxt(eigenval_file)

# Calculate variance explained
total_var = np.sum(eigenvalues)
var_explained = (eigenvalues / total_var) * 100
cumulative_var = np.cumsum(var_explained)

# Create scree plot
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))

# Bar plot of variance explained
pcs = range(1, len(eigenvalues) + 1)
ax1.bar(pcs, var_explained, color='steelblue', edgecolor='black')
ax1.set_xlabel('Principal Component')
ax1.set_ylabel('Variance Explained (%)')
ax1.set_title('PCA Scree Plot')
ax1.set_xticks(pcs)

# Cumulative variance line plot
ax2.plot(pcs, cumulative_var, 'o-', color='darkred', linewidth=2, markersize=6)
ax2.axhline(y=80, color='gray', linestyle='--', label='80% threshold')
ax2.set_xlabel('Principal Component')
ax2.set_ylabel('Cumulative Variance Explained (%)')
ax2.set_title('Cumulative Variance')
ax2.set_xticks(pcs)
ax2.legend()

plt.tight_layout()
plt.savefig(os.path.join(results_dir, 'scree_plot.png'), dpi=150, bbox_inches='tight')
plt.close()

# Save variance summary
with open(os.path.join(results_dir, 'variance_summary.txt'), 'w') as f:
    f.write("PC\tVariance_Explained(%)\tCumulative(%)\n")
    for i, (ve, cv) in enumerate(zip(var_explained, cumulative_var)):
        f.write(f"PC{i+1}\t{ve:.2f}\t{cv:.2f}\n")

print(f"Scree plot saved to {os.path.join(results_dir, 'scree_plot.png')}")
EOF
    
    echo "Scree plot generated" | tee -a "$LOG_FILE"
}

#-------------------------------------------------------------------------------
# Main execution
#-------------------------------------------------------------------------------
main() {
    echo "" | tee -a "$LOG_FILE"
    
    # Check for required tools
    for tool in plink python3; do
        if ! command -v "$tool" &> /dev/null; then
            echo "ERROR: Required tool '$tool' not found." | tee -a "$LOG_FILE"
            exit 1
        fi
    done
    
    # Check inputs
    check_inputs
    
    # QC filtering
    qc_prefix=$(qc_filter)
    
    # LD pruning
    pruned_prefix=$(ld_prune "$qc_prefix")
    
    # Run PCA
    run_pca "$pruned_prefix"
    
    # Generate plots
    generate_scree_plot
    
    echo "" | tee -a "$LOG_FILE"
    echo "============================================" | tee -a "$LOG_FILE"
    echo "PCA analysis completed successfully!" | tee -a "$LOG_FILE"
    echo "Results in: ${RESULTS_DIR}" | tee -a "$LOG_FILE"
    echo "Finished: $(date)" | tee -a "$LOG_FILE"
    echo "============================================" | tee -a "$LOG_FILE"
}

# Run main function
main "$@"
