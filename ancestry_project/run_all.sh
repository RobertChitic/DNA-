#!/bin/bash
#===============================================================================
# Master Script: run_all.sh
# Description: Execute complete ancestry analysis pipeline
#
# Usage: ./run_all.sh <ancestry_file.txt> [user_id]
#
# This script runs all analysis steps in sequence:
# 1. Download reference data (if not present)
# 2. Convert AncestryDNA file to PLINK format
# 3. Merge user data with references
# 4. Run PCA analysis
# 5. Run ADMIXTURE analysis
# 6. Calculate ancestry probabilities
# 7. Generate ensemble results
# 8. Create visualizations and HTML report
#
# Requirements:
# - Conda environment with all dependencies
# - 75GB+ disk space (for reference data)
# - 4-8 hours runtime (depending on system)
#===============================================================================

set -euo pipefail

# Configuration
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_DIR="$(dirname "$SCRIPT_DIR")"

# Parameters
ANCESTRY_FILE="${1:-}"
USER_ID="${2:-USER001}"

# Output directories
DATA_DIR="${PROJECT_DIR}/data"
RAW_DIR="${DATA_DIR}/raw"
PROCESSED_DIR="${DATA_DIR}/processed"
REF_DIR="${PROJECT_DIR}/reference"
RESULTS_DIR="${PROJECT_DIR}/results"

# Default parameters (can be overridden with environment variables)
N_PCS="${N_PCS:-10}"
BANDWIDTH="${BANDWIDTH:-0.5}"
MIN_K="${MIN_K:-2}"
MAX_K="${MAX_K:-15}"
BOOTSTRAP="${BOOTSTRAP:-0}"
PCA_WEIGHT="${PCA_WEIGHT:-0.3}"
ADMIXTURE_WEIGHT="${ADMIXTURE_WEIGHT:-0.7}"

# Log file
LOG_FILE="${RESULTS_DIR}/pipeline.log"

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

#-------------------------------------------------------------------------------
# Utility Functions
#-------------------------------------------------------------------------------

log_info() {
    echo -e "${BLUE}[INFO]${NC} $(date '+%Y-%m-%d %H:%M:%S') - $1" | tee -a "$LOG_FILE"
}

log_success() {
    echo -e "${GREEN}[SUCCESS]${NC} $(date '+%Y-%m-%d %H:%M:%S') - $1" | tee -a "$LOG_FILE"
}

log_warning() {
    echo -e "${YELLOW}[WARNING]${NC} $(date '+%Y-%m-%d %H:%M:%S') - $1" | tee -a "$LOG_FILE"
}

log_error() {
    echo -e "${RED}[ERROR]${NC} $(date '+%Y-%m-%d %H:%M:%S') - $1" | tee -a "$LOG_FILE"
}

check_step_completed() {
    local marker_file="$1"
    if [[ -f "$marker_file" ]]; then
        return 0
    fi
    return 1
}

mark_step_completed() {
    local marker_file="$1"
    touch "$marker_file"
}

#-------------------------------------------------------------------------------
# Print Banner
#-------------------------------------------------------------------------------
print_banner() {
    echo ""
    echo "=============================================="
    echo "    🧬 DIY Ancestry Analysis Pipeline 🧬    "
    echo "=============================================="
    echo ""
    echo "Project Directory: $PROJECT_DIR"
    echo "Ancestry File: $ANCESTRY_FILE"
    echo "User ID: $USER_ID"
    echo ""
    echo "Parameters:"
    echo "  N_PCS: $N_PCS"
    echo "  BANDWIDTH: $BANDWIDTH"
    echo "  K Range: $MIN_K - $MAX_K"
    echo "  Bootstrap: $BOOTSTRAP"
    echo "  Weights: PCA=$PCA_WEIGHT, ADMIXTURE=$ADMIXTURE_WEIGHT"
    echo ""
    echo "=============================================="
    echo ""
}

#-------------------------------------------------------------------------------
# Check Prerequisites
#-------------------------------------------------------------------------------
check_prerequisites() {
    log_info "Checking prerequisites..."
    
    # Check for required tools
    local required_tools=("plink" "python3" "wget")
    local missing_tools=()
    
    for tool in "${required_tools[@]}"; do
        if ! command -v "$tool" &> /dev/null; then
            missing_tools+=("$tool")
        fi
    done
    
    if [[ ${#missing_tools[@]} -gt 0 ]]; then
        log_error "Missing required tools: ${missing_tools[*]}"
        log_info "Please install missing tools or activate conda environment:"
        log_info "  conda env create -f ${PROJECT_DIR}/environment.yml"
        log_info "  conda activate ancestry_analysis"
        exit 1
    fi
    
    # Check Python packages
    python3 -c "import pandas; import numpy; import scipy; import matplotlib" 2>/dev/null || {
        log_error "Missing Python packages. Please install with:"
        log_info "  pip install pandas numpy scipy matplotlib seaborn plotly"
        exit 1
    }
    
    log_success "All prerequisites satisfied"
}

#-------------------------------------------------------------------------------
# Step 1: Download Reference Data
#-------------------------------------------------------------------------------
step_download_data() {
    local marker="${RESULTS_DIR}/.step1_complete"
    
    if check_step_completed "$marker"; then
        log_info "Step 1: Reference data already downloaded, skipping"
        return 0
    fi
    
    log_info "Step 1: Downloading reference data..."
    
    if [[ -f "${RAW_DIR}/1kg_euro_merged.vcf.gz" ]] && [[ -f "${RAW_DIR}/hgdp_euro.vcf.gz" ]]; then
        log_info "Reference VCF files already exist, skipping download"
    else
        bash "${SCRIPT_DIR}/01_download_data.sh" 2>&1 | tee -a "$LOG_FILE"
    fi
    
    mark_step_completed "$marker"
    log_success "Step 1 complete"
}

#-------------------------------------------------------------------------------
# Step 2: Convert AncestryDNA File
#-------------------------------------------------------------------------------
step_convert_ancestry() {
    local marker="${RESULTS_DIR}/.step2_complete"
    
    if check_step_completed "$marker" && [[ -f "${PROCESSED_DIR}/user.bed" ]]; then
        log_info "Step 2: User data already converted, skipping"
        return 0
    fi
    
    log_info "Step 2: Converting AncestryDNA file to PLINK format..."
    
    if [[ -z "$ANCESTRY_FILE" ]] || [[ ! -f "$ANCESTRY_FILE" ]]; then
        log_error "Ancestry file not provided or not found: $ANCESTRY_FILE"
        log_info "Usage: ./run_all.sh <ancestry_file.txt> [user_id]"
        exit 1
    fi
    
    python3 "${SCRIPT_DIR}/02_convert_ancestry_file.py" \
        --input "$ANCESTRY_FILE" \
        --output "${PROCESSED_DIR}/user" \
        --sample-id "$USER_ID" \
        2>&1 | tee -a "$LOG_FILE"
    
    mark_step_completed "$marker"
    log_success "Step 2 complete"
}

#-------------------------------------------------------------------------------
# Step 3: Merge Datasets
#-------------------------------------------------------------------------------
step_merge_data() {
    local marker="${RESULTS_DIR}/.step3_complete"
    
    if check_step_completed "$marker" && [[ -f "${PROCESSED_DIR}/merged_all.bed" ]]; then
        log_info "Step 3: Data already merged, skipping"
        return 0
    fi
    
    log_info "Step 3: Merging user data with reference panels..."
    
    bash "${SCRIPT_DIR}/03_merge_references.sh" "${PROCESSED_DIR}/user" \
        2>&1 | tee -a "$LOG_FILE"
    
    mark_step_completed "$marker"
    log_success "Step 3 complete"
}

#-------------------------------------------------------------------------------
# Step 4: Run PCA
#-------------------------------------------------------------------------------
step_run_pca() {
    local marker="${RESULTS_DIR}/.step4_complete"
    
    if check_step_completed "$marker" && [[ -f "${RESULTS_DIR}/pca/eigenvec.txt" ]]; then
        log_info "Step 4: PCA already completed, skipping"
        return 0
    fi
    
    log_info "Step 4: Running PCA analysis..."
    
    bash "${SCRIPT_DIR}/04_run_pca.sh" "${PROCESSED_DIR}/merged_all" "$N_PCS" \
        2>&1 | tee -a "$LOG_FILE"
    
    mark_step_completed "$marker"
    log_success "Step 4 complete"
}

#-------------------------------------------------------------------------------
# Step 5: Run ADMIXTURE
#-------------------------------------------------------------------------------
step_run_admixture() {
    local marker="${RESULTS_DIR}/.step5_complete"
    
    if check_step_completed "$marker" && [[ -f "${RESULTS_DIR}/admixture/optimal_k.txt" ]]; then
        log_info "Step 5: ADMIXTURE already completed, skipping"
        return 0
    fi
    
    log_info "Step 5: Running ADMIXTURE analysis..."
    
    bash "${SCRIPT_DIR}/05_run_admixture.sh" \
        "${RESULTS_DIR}/pca/ld_pruned" "$MIN_K" "$MAX_K" \
        2>&1 | tee -a "$LOG_FILE"
    
    mark_step_completed "$marker"
    log_success "Step 5 complete"
}

#-------------------------------------------------------------------------------
# Step 6: Calculate PCA Probabilities
#-------------------------------------------------------------------------------
step_pca_probabilities() {
    local marker="${RESULTS_DIR}/.step6_complete"
    local output_file="${RESULTS_DIR}/probabilities/pca_probabilities.csv"
    
    if check_step_completed "$marker" && [[ -f "$output_file" ]]; then
        log_info "Step 6: PCA probabilities already calculated, skipping"
        return 0
    fi
    
    log_info "Step 6: Calculating PCA-based ancestry probabilities..."
    
    mkdir -p "${RESULTS_DIR}/probabilities"
    
    python3 "${SCRIPT_DIR}/06_pca_probabilities.py" \
        --eigenvec "${RESULTS_DIR}/pca/eigenvec.txt" \
        --metadata "${REF_DIR}/sample_metadata.tsv" \
        --output "$output_file" \
        --n-pcs "$N_PCS" \
        --bandwidth "$BANDWIDTH" \
        --bootstrap "$BOOTSTRAP" \
        --user-id "$USER_ID" \
        2>&1 | tee -a "$LOG_FILE"
    
    mark_step_completed "$marker"
    log_success "Step 6 complete"
}

#-------------------------------------------------------------------------------
# Step 7: Calculate ADMIXTURE Probabilities
#-------------------------------------------------------------------------------
step_admixture_probabilities() {
    local marker="${RESULTS_DIR}/.step7_complete"
    local output_file="${RESULTS_DIR}/probabilities/admixture_probabilities.csv"
    
    if check_step_completed "$marker" && [[ -f "$output_file" ]]; then
        log_info "Step 7: ADMIXTURE probabilities already calculated, skipping"
        return 0
    fi
    
    log_info "Step 7: Calculating ADMIXTURE-based ancestry probabilities..."
    
    # Get optimal K
    local optimal_k=$(cat "${RESULTS_DIR}/admixture/optimal_k.txt")
    local q_file="${RESULTS_DIR}/admixture/ld_pruned.${optimal_k}.Q"
    local fam_file="${RESULTS_DIR}/admixture/samples.fam"
    
    python3 "${SCRIPT_DIR}/07_admixture_probabilities.py" \
        --q-file "$q_file" \
        --fam-file "$fam_file" \
        --metadata "${REF_DIR}/sample_metadata.tsv" \
        --output "$output_file" \
        --bootstrap "$BOOTSTRAP" \
        --user-id "$USER_ID" \
        2>&1 | tee -a "$LOG_FILE"
    
    mark_step_completed "$marker"
    log_success "Step 7 complete"
}

#-------------------------------------------------------------------------------
# Step 8: Generate Ensemble Results
#-------------------------------------------------------------------------------
step_ensemble() {
    local marker="${RESULTS_DIR}/.step8_complete"
    local output_file="${RESULTS_DIR}/FINAL_ensemble.csv"
    
    if check_step_completed "$marker" && [[ -f "$output_file" ]]; then
        log_info "Step 8: Ensemble results already generated, skipping"
        return 0
    fi
    
    log_info "Step 8: Generating ensemble ancestry estimates..."
    
    python3 "${SCRIPT_DIR}/09_ensemble.py" \
        --pca-probs "${RESULTS_DIR}/probabilities/pca_probabilities.csv" \
        --admixture-probs "${RESULTS_DIR}/probabilities/admixture_probabilities.csv" \
        --output "$output_file" \
        --pca-weight "$PCA_WEIGHT" \
        --admixture-weight "$ADMIXTURE_WEIGHT" \
        2>&1 | tee -a "$LOG_FILE"
    
    mark_step_completed "$marker"
    log_success "Step 8 complete"
}

#-------------------------------------------------------------------------------
# Step 9: Generate Visualizations
#-------------------------------------------------------------------------------
step_visualizations() {
    local marker="${RESULTS_DIR}/.step9_complete"
    
    if check_step_completed "$marker"; then
        log_info "Step 9: Visualizations already generated, skipping"
        return 0
    fi
    
    log_info "Step 9: Generating visualizations..."
    
    python3 "${SCRIPT_DIR}/08_visualizations.py" \
        --results-dir "$RESULTS_DIR" \
        --user-id "$USER_ID" \
        2>&1 | tee -a "$LOG_FILE"
    
    mark_step_completed "$marker"
    log_success "Step 9 complete"
}

#-------------------------------------------------------------------------------
# Step 10: Generate HTML Report
#-------------------------------------------------------------------------------
step_generate_report() {
    local marker="${RESULTS_DIR}/.step10_complete"
    local output_file="${RESULTS_DIR}/interactive_dashboard.html"
    
    if check_step_completed "$marker" && [[ -f "$output_file" ]]; then
        log_info "Step 10: HTML report already generated, skipping"
        return 0
    fi
    
    log_info "Step 10: Generating interactive HTML report..."
    
    python3 "${SCRIPT_DIR}/10_generate_report.py" \
        --results-dir "$RESULTS_DIR" \
        --output "$output_file" \
        2>&1 | tee -a "$LOG_FILE"
    
    mark_step_completed "$marker"
    log_success "Step 10 complete"
}

#-------------------------------------------------------------------------------
# Print Final Summary
#-------------------------------------------------------------------------------
print_summary() {
    echo ""
    echo "=============================================="
    echo "        🎉 Analysis Complete! 🎉            "
    echo "=============================================="
    echo ""
    echo "Results saved to: ${RESULTS_DIR}"
    echo ""
    echo "Output files:"
    echo "  📊 ${RESULTS_DIR}/FINAL_ensemble.csv"
    echo "  📈 ${RESULTS_DIR}/pca_plot.png"
    echo "  📈 ${RESULTS_DIR}/admixture_barplot.png"
    echo "  📈 ${RESULTS_DIR}/comparison_plot.png"
    echo "  🌐 ${RESULTS_DIR}/interactive_dashboard.html"
    echo ""
    
    # Print top results if available
    if [[ -f "${RESULTS_DIR}/FINAL_ensemble.csv" ]]; then
        echo "Top Ancestry Results:"
        echo "---------------------"
        head -11 "${RESULTS_DIR}/FINAL_ensemble.csv" | column -t -s,
    fi
    
    echo ""
    echo "=============================================="
    echo "To view detailed results, open:"
    echo "  ${RESULTS_DIR}/interactive_dashboard.html"
    echo "=============================================="
}

#-------------------------------------------------------------------------------
# Clean up function (for reruns)
#-------------------------------------------------------------------------------
clean_markers() {
    log_info "Cleaning step markers for fresh run..."
    rm -f "${RESULTS_DIR}"/.step*_complete
}

#-------------------------------------------------------------------------------
# Main Execution
#-------------------------------------------------------------------------------
main() {
    # Create necessary directories
    mkdir -p "$RAW_DIR" "$PROCESSED_DIR" "$REF_DIR" "$RESULTS_DIR"
    
    # Initialize log
    echo "Pipeline started: $(date)" > "$LOG_FILE"
    
    # Print banner
    print_banner
    
    # Check for --clean flag
    if [[ "${1:-}" == "--clean" ]]; then
        clean_markers
        shift
        ANCESTRY_FILE="${1:-}"
        USER_ID="${2:-USER001}"
    fi
    
    # Check prerequisites
    check_prerequisites
    
    # Execute pipeline steps
    step_download_data
    step_convert_ancestry
    step_merge_data
    step_run_pca
    step_run_admixture
    step_pca_probabilities
    step_admixture_probabilities
    step_ensemble
    step_visualizations
    step_generate_report
    
    # Print summary
    print_summary
    
    log_info "Pipeline completed successfully!"
}

# Show usage if no arguments
if [[ $# -eq 0 ]] && [[ ! -f "${RAW_DIR}/1kg_euro_merged.vcf.gz" ]]; then
    echo "Usage: $0 <ancestry_file.txt> [user_id]"
    echo ""
    echo "Arguments:"
    echo "  ancestry_file.txt  - Your AncestryDNA raw data file"
    echo "  user_id           - Identifier for your sample (default: USER001)"
    echo ""
    echo "Options:"
    echo "  --clean           - Remove step markers and rerun from scratch"
    echo ""
    echo "Environment variables (optional):"
    echo "  N_PCS=10          - Number of principal components"
    echo "  BANDWIDTH=0.5     - Gaussian kernel bandwidth"
    echo "  MIN_K=2           - Minimum K for ADMIXTURE"
    echo "  MAX_K=15          - Maximum K for ADMIXTURE"
    echo "  BOOTSTRAP=0       - Bootstrap iterations (0 to disable)"
    echo ""
    echo "Example:"
    echo "  ./run_all.sh ~/Downloads/AncestryDNA.txt MyName"
    echo ""
    exit 0
fi

# Run main
main "$@"
