#!/bin/bash
#===============================================================================
# Script: 05_run_admixture.sh
# Description: Run ADMIXTURE analysis with cross-validation
#
# Usage: ./05_run_admixture.sh [input_prefix] [min_k] [max_k]
#
# Input:
#   - LD-pruned PLINK files from 04_run_pca.sh
#
# Output:
#   - results/admixture/*.P files (population allele frequencies)
#   - results/admixture/*.Q files (individual ancestry proportions)
#   - results/admixture/cv_errors.txt (cross-validation errors)
#   - results/admixture/cv_plot.png (CV error plot)
#===============================================================================

set -euo pipefail

# Configuration
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_DIR="$(dirname "$SCRIPT_DIR")"
PROCESSED_DIR="${PROJECT_DIR}/data/processed"
PCA_DIR="${PROJECT_DIR}/results/pca"
RESULTS_DIR="${PROJECT_DIR}/results/admixture"
LOG_FILE="${PROJECT_DIR}/results/05_admixture.log"

# Parameters
INPUT_PREFIX="${1:-${PCA_DIR}/ld_pruned}"
MIN_K="${2:-2}"
MAX_K="${3:-15}"
N_THREADS="${N_THREADS:-4}"

# Create directories
mkdir -p "$RESULTS_DIR" "$(dirname "$LOG_FILE")"

# Initialize log
echo "============================================" | tee "$LOG_FILE"
echo "Ancestry Analysis - ADMIXTURE" | tee -a "$LOG_FILE"
echo "Started: $(date)" | tee -a "$LOG_FILE"
echo "============================================" | tee -a "$LOG_FILE"
echo "Input: ${INPUT_PREFIX}" | tee -a "$LOG_FILE"
echo "K range: ${MIN_K} to ${MAX_K}" | tee -a "$LOG_FILE"
echo "Threads: ${N_THREADS}" | tee -a "$LOG_FILE"

#-------------------------------------------------------------------------------
# Check input files
#-------------------------------------------------------------------------------
check_inputs() {
    if [[ ! -f "${INPUT_PREFIX}.bed" ]]; then
        echo "ERROR: Input PLINK files not found: ${INPUT_PREFIX}" | tee -a "$LOG_FILE"
        echo "Please run 04_run_pca.sh first." | tee -a "$LOG_FILE"
        exit 1
    fi
    
    local n_variants=$(wc -l < "${INPUT_PREFIX}.bim")
    local n_samples=$(wc -l < "${INPUT_PREFIX}.fam")
    echo "Input dataset: ${n_variants} variants, ${n_samples} samples" | tee -a "$LOG_FILE"
}

#-------------------------------------------------------------------------------
# Check for ADMIXTURE
#-------------------------------------------------------------------------------
check_admixture() {
    if ! command -v admixture &> /dev/null; then
        echo "ADMIXTURE not found. Attempting to install..." | tee -a "$LOG_FILE"
        
        # Try conda installation
        if command -v conda &> /dev/null; then
            conda install -c bioconda admixture -y >> "$LOG_FILE" 2>&1 || true
        fi
        
        # Check again
        if ! command -v admixture &> /dev/null; then
            echo "ERROR: ADMIXTURE not found. Please install ADMIXTURE 1.3.0:" | tee -a "$LOG_FILE"
            echo "  wget http://dalexander.github.io/admixture/binaries/admixture_linux-1.3.0.tar.gz" | tee -a "$LOG_FILE"
            echo "  tar xzf admixture_linux-1.3.0.tar.gz" | tee -a "$LOG_FILE"
            echo "  sudo mv admixture /usr/local/bin/" | tee -a "$LOG_FILE"
            exit 1
        fi
    fi
    
    echo "ADMIXTURE version: $(admixture 2>&1 | head -1)" | tee -a "$LOG_FILE"
}

#-------------------------------------------------------------------------------
# Run ADMIXTURE for all K values
#-------------------------------------------------------------------------------
run_admixture() {
    echo "" | tee -a "$LOG_FILE"
    echo "=== Running ADMIXTURE ===" | tee -a "$LOG_FILE"
    
    # Get basename for input file
    local input_basename=$(basename "$INPUT_PREFIX")
    
    # Create working directory
    local work_dir="${RESULTS_DIR}/work"
    mkdir -p "$work_dir"
    
    # Copy input files to working directory (ADMIXTURE writes to current directory)
    cp "${INPUT_PREFIX}.bed" "${work_dir}/${input_basename}.bed"
    cp "${INPUT_PREFIX}.bim" "${work_dir}/${input_basename}.bim"
    cp "${INPUT_PREFIX}.fam" "${work_dir}/${input_basename}.fam"
    
    # Change to working directory
    cd "$work_dir"
    
    # Initialize CV error file
    local cv_file="${RESULTS_DIR}/cv_errors.txt"
    echo "K	CV_Error" > "$cv_file"
    
    # Run ADMIXTURE for each K
    for k in $(seq $MIN_K $MAX_K); do
        echo "Running K=${k}..." | tee -a "$LOG_FILE"
        
        local q_file="${input_basename}.${k}.Q"
        local p_file="${input_basename}.${k}.P"
        
        if [[ ! -f "$q_file" ]]; then
            # Run ADMIXTURE with cross-validation
            admixture --cv -j${N_THREADS} "${input_basename}.bed" $k 2>&1 | tee -a "$LOG_FILE"
            
            # Extract CV error
            local cv_error=$(grep "CV error" "${work_dir}/admixture_log_k${k}.txt" 2>/dev/null | awk '{print $NF}' || grep -o "CV error.*" *.log 2>/dev/null | tail -1 | awk '{print $NF}')
            
            if [[ -z "$cv_error" ]]; then
                # Try to get CV error from stdout saved to log
                cv_error=$(grep "CV error (K=${k})" "$LOG_FILE" | tail -1 | awk '{print $NF}')
            fi
            
            if [[ -n "$cv_error" ]]; then
                echo "${k}	${cv_error}" >> "$cv_file"
            fi
        else
            echo "K=${k} results already exist, skipping." | tee -a "$LOG_FILE"
        fi
        
        # Copy results to main results directory
        [[ -f "$q_file" ]] && cp "$q_file" "$RESULTS_DIR/"
        [[ -f "$p_file" ]] && cp "$p_file" "$RESULTS_DIR/"
    done
    
    # Return to original directory
    cd - > /dev/null
    
    echo "ADMIXTURE runs completed" | tee -a "$LOG_FILE"
}

#-------------------------------------------------------------------------------
# Find optimal K
#-------------------------------------------------------------------------------
find_optimal_k() {
    echo "" | tee -a "$LOG_FILE"
    echo "=== Finding Optimal K ===" | tee -a "$LOG_FILE"
    
    local cv_file="${RESULTS_DIR}/cv_errors.txt"
    
    if [[ ! -f "$cv_file" ]] || [[ $(wc -l < "$cv_file") -le 1 ]]; then
        echo "WARNING: CV errors not available, using K=${MIN_K}" | tee -a "$LOG_FILE"
        echo "$MIN_K" > "${RESULTS_DIR}/optimal_k.txt"
        return
    fi
    
    # Find K with minimum CV error
    local optimal_k=$(tail -n +2 "$cv_file" | sort -t$'\t' -k2 -n | head -1 | cut -f1)
    
    echo "Optimal K: ${optimal_k}" | tee -a "$LOG_FILE"
    echo "$optimal_k" > "${RESULTS_DIR}/optimal_k.txt"
}

#-------------------------------------------------------------------------------
# Generate CV error plot
#-------------------------------------------------------------------------------
generate_cv_plot() {
    echo "" | tee -a "$LOG_FILE"
    echo "=== Generating CV Error Plot ===" | tee -a "$LOG_FILE"
    
    python3 << 'EOF'
import matplotlib.pyplot as plt
import pandas as pd
import os

# Read CV errors
results_dir = os.environ.get('RESULTS_DIR', 'results/admixture')
cv_file = os.path.join(results_dir, 'cv_errors.txt')

if not os.path.exists(cv_file):
    print(f"CV file not found: {cv_file}")
    exit(0)

df = pd.read_csv(cv_file, sep='\t')

if df.empty or len(df) < 2:
    print("Not enough CV data to plot")
    exit(0)

# Find optimal K
optimal_k = df.loc[df['CV_Error'].idxmin(), 'K']

# Create plot
plt.figure(figsize=(10, 6))
plt.plot(df['K'], df['CV_Error'], 'o-', color='steelblue', linewidth=2, markersize=8)
plt.axvline(x=optimal_k, color='red', linestyle='--', label=f'Optimal K={int(optimal_k)}')
plt.scatter([optimal_k], [df[df['K']==optimal_k]['CV_Error'].values[0]], 
            color='red', s=150, zorder=5, marker='*')

plt.xlabel('K (Number of Ancestral Populations)', fontsize=12)
plt.ylabel('Cross-Validation Error', fontsize=12)
plt.title('ADMIXTURE Cross-Validation Error', fontsize=14)
plt.legend()
plt.xticks(df['K'])
plt.grid(True, alpha=0.3)

plt.tight_layout()
plt.savefig(os.path.join(results_dir, 'cv_plot.png'), dpi=150, bbox_inches='tight')
plt.close()

print(f"CV plot saved with optimal K={int(optimal_k)}")
EOF
    
    echo "CV error plot generated" | tee -a "$LOG_FILE"
}

#-------------------------------------------------------------------------------
# Create sample-population mapping for interpretation
#-------------------------------------------------------------------------------
create_pop_mapping() {
    echo "" | tee -a "$LOG_FILE"
    echo "=== Creating Population Mapping ===" | tee -a "$LOG_FILE"
    
    # Copy FAM file sample info
    local fam_file="${INPUT_PREFIX}.fam"
    local pop_file="${INPUT_PREFIX}.pop"
    
    if [[ -f "$pop_file" ]]; then
        cp "$pop_file" "${RESULTS_DIR}/sample_populations.txt"
    else
        # Create from FAM file
        awk '{print $1}' "$fam_file" > "${RESULTS_DIR}/sample_ids.txt"
    fi
    
    # Copy FAM for reference
    cp "$fam_file" "${RESULTS_DIR}/samples.fam"
}

#-------------------------------------------------------------------------------
# Main execution
#-------------------------------------------------------------------------------
main() {
    echo "" | tee -a "$LOG_FILE"
    
    # Export RESULTS_DIR for Python script
    export RESULTS_DIR
    
    # Check for required tools
    if ! command -v python3 &> /dev/null; then
        echo "ERROR: Python3 not found." | tee -a "$LOG_FILE"
        exit 1
    fi
    
    # Check inputs
    check_inputs
    
    # Check ADMIXTURE installation
    check_admixture
    
    # Run ADMIXTURE
    run_admixture
    
    # Find optimal K
    find_optimal_k
    
    # Generate plots
    generate_cv_plot
    
    # Create mappings
    create_pop_mapping
    
    echo "" | tee -a "$LOG_FILE"
    echo "============================================" | tee -a "$LOG_FILE"
    echo "ADMIXTURE analysis completed successfully!" | tee -a "$LOG_FILE"
    echo "Results in: ${RESULTS_DIR}" | tee -a "$LOG_FILE"
    echo "Optimal K: $(cat ${RESULTS_DIR}/optimal_k.txt)" | tee -a "$LOG_FILE"
    echo "Finished: $(date)" | tee -a "$LOG_FILE"
    echo "============================================" | tee -a "$LOG_FILE"
}

# Run main function
main "$@"
