#!/bin/bash
#===============================================================================
# Script: 03_merge_references.sh
# Description: Merge reference datasets with user data, handling strand flips
#
# Usage: ./03_merge_references.sh <user_prefix>
#
# Input:
#   - User PLINK files (from 02_convert_ancestry_file.py)
#   - Reference PLINK files (from 01_download_data.sh)
#
# Output:
#   - data/processed/merged_all.bed/bim/fam (merged dataset)
#===============================================================================

set -euo pipefail

# Configuration
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_DIR="$(dirname "$SCRIPT_DIR")"
RAW_DIR="${PROJECT_DIR}/data/raw"
PROCESSED_DIR="${PROJECT_DIR}/data/processed"
REF_DIR="${PROJECT_DIR}/reference"
LOG_FILE="${PROJECT_DIR}/results/03_merge.log"

# Default user prefix
USER_PREFIX="${1:-${PROCESSED_DIR}/user}"

# Create directories
mkdir -p "$PROCESSED_DIR" "$(dirname "$LOG_FILE")"

# Initialize log
echo "============================================" | tee "$LOG_FILE"
echo "Ancestry Analysis - Data Merge" | tee -a "$LOG_FILE"
echo "Started: $(date)" | tee -a "$LOG_FILE"
echo "============================================" | tee -a "$LOG_FILE"

#-------------------------------------------------------------------------------
# Function: check_plink_files
# Description: Verify PLINK binary files exist
#-------------------------------------------------------------------------------
check_plink_files() {
    local prefix="$1"
    local name="$2"
    
    if [[ ! -f "${prefix}.bed" ]] || [[ ! -f "${prefix}.bim" ]] || [[ ! -f "${prefix}.fam" ]]; then
        echo "ERROR: Missing PLINK files for ${name}: ${prefix}" | tee -a "$LOG_FILE"
        return 1
    fi
    
    local n_variants=$(wc -l < "${prefix}.bim")
    local n_samples=$(wc -l < "${prefix}.fam")
    echo "${name}: ${n_variants} variants, ${n_samples} samples" | tee -a "$LOG_FILE"
    
    return 0
}

#-------------------------------------------------------------------------------
# Function: convert_vcf_to_plink
# Description: Convert VCF to PLINK binary format
#-------------------------------------------------------------------------------
convert_vcf_to_plink() {
    local vcf_file="$1"
    local output_prefix="$2"
    
    if [[ ! -f "${output_prefix}.bed" ]]; then
        echo "Converting ${vcf_file} to PLINK format..." | tee -a "$LOG_FILE"
        
        plink --vcf "$vcf_file" \
            --make-bed \
            --out "$output_prefix" \
            --double-id \
            --allow-extra-chr \
            --snps-only just-acgt \
            --const-fid 0 \
            >> "$LOG_FILE" 2>&1
    else
        echo "PLINK files for ${output_prefix} already exist, skipping conversion." | tee -a "$LOG_FILE"
    fi
}

#-------------------------------------------------------------------------------
# Function: handle_strand_flips
# Description: Identify and handle strand flips between datasets
#-------------------------------------------------------------------------------
handle_strand_flips() {
    local ref_prefix="$1"
    local target_prefix="$2"
    local output_prefix="$3"
    
    echo "Checking for strand flips..." | tee -a "$LOG_FILE"
    
    # Create temp directory for intermediate files
    local temp_dir=$(mktemp -d)
    
    # Find common SNPs
    cut -f2 "${ref_prefix}.bim" | sort > "${temp_dir}/ref_snps.txt"
    cut -f2 "${target_prefix}.bim" | sort > "${temp_dir}/target_snps.txt"
    comm -12 "${temp_dir}/ref_snps.txt" "${temp_dir}/target_snps.txt" > "${temp_dir}/common_snps.txt"
    
    local n_common=$(wc -l < "${temp_dir}/common_snps.txt")
    echo "Found ${n_common} common SNPs" | tee -a "$LOG_FILE"
    
    # Extract common SNPs from both datasets
    plink --bfile "$ref_prefix" \
        --extract "${temp_dir}/common_snps.txt" \
        --make-bed \
        --out "${temp_dir}/ref_common" \
        >> "$LOG_FILE" 2>&1
    
    plink --bfile "$target_prefix" \
        --extract "${temp_dir}/common_snps.txt" \
        --make-bed \
        --out "${temp_dir}/target_common" \
        >> "$LOG_FILE" 2>&1
    
    # Identify strand flips by comparing alleles
    # Create allele comparison file
    awk '{print $2, $5, $6}' "${temp_dir}/ref_common.bim" | sort -k1 > "${temp_dir}/ref_alleles.txt"
    awk '{print $2, $5, $6}' "${temp_dir}/target_common.bim" | sort -k1 > "${temp_dir}/target_alleles.txt"
    
    # Find SNPs that need flipping (complementary alleles)
    python3 << EOF > "${temp_dir}/flip_snps.txt"
import sys

complement = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}

ref_alleles = {}
with open('${temp_dir}/ref_alleles.txt') as f:
    for line in f:
        parts = line.strip().split()
        if len(parts) >= 3:
            ref_alleles[parts[0]] = (parts[1], parts[2])

flip_snps = []
with open('${temp_dir}/target_alleles.txt') as f:
    for line in f:
        parts = line.strip().split()
        if len(parts) >= 3 and parts[0] in ref_alleles:
            snp = parts[0]
            t_a1, t_a2 = parts[1], parts[2]
            r_a1, r_a2 = ref_alleles[snp]
            
            # Check if strand flip needed
            t_a1_comp = complement.get(t_a1, t_a1)
            t_a2_comp = complement.get(t_a2, t_a2)
            
            # If complementary alleles match reference
            if {t_a1_comp, t_a2_comp} == {r_a1, r_a2}:
                print(snp)
                
EOF

    local n_flip=$(wc -l < "${temp_dir}/flip_snps.txt")
    echo "Found ${n_flip} SNPs requiring strand flip" | tee -a "$LOG_FILE"
    
    # Flip strands if needed
    if [[ $n_flip -gt 0 ]]; then
        plink --bfile "${temp_dir}/target_common" \
            --flip "${temp_dir}/flip_snps.txt" \
            --make-bed \
            --out "${temp_dir}/target_flipped" \
            >> "$LOG_FILE" 2>&1
        
        cp "${temp_dir}/target_flipped.bed" "${output_prefix}.bed"
        cp "${temp_dir}/target_flipped.bim" "${output_prefix}.bim"
        cp "${temp_dir}/target_flipped.fam" "${output_prefix}.fam"
    else
        cp "${temp_dir}/target_common.bed" "${output_prefix}.bed"
        cp "${temp_dir}/target_common.bim" "${output_prefix}.bim"
        cp "${temp_dir}/target_common.fam" "${output_prefix}.fam"
    fi
    
    # Clean up
    rm -rf "$temp_dir"
    
    echo "Strand flip handling complete" | tee -a "$LOG_FILE"
}

#-------------------------------------------------------------------------------
# Function: merge_datasets
# Description: Merge multiple PLINK datasets
#-------------------------------------------------------------------------------
merge_datasets() {
    local output_prefix="$1"
    shift
    local datasets=("$@")
    
    echo "Merging ${#datasets[@]} datasets..." | tee -a "$LOG_FILE"
    
    if [[ ${#datasets[@]} -lt 2 ]]; then
        echo "ERROR: Need at least 2 datasets to merge" | tee -a "$LOG_FILE"
        return 1
    fi
    
    # Create merge list (all datasets except first)
    local merge_list="${PROCESSED_DIR}/merge_list.txt"
    rm -f "$merge_list"
    
    for ((i=1; i<${#datasets[@]}; i++)); do
        echo "${datasets[$i]}" >> "$merge_list"
    done
    
    # Perform merge
    plink --bfile "${datasets[0]}" \
        --merge-list "$merge_list" \
        --make-bed \
        --out "$output_prefix" \
        --allow-no-sex \
        >> "$LOG_FILE" 2>&1
    
    # Check if merge failed due to strand issues
    if [[ -f "${output_prefix}-merge.missnp" ]]; then
        echo "Handling problematic SNPs from merge..." | tee -a "$LOG_FILE"
        
        # Remove problematic SNPs from all datasets and retry
        local temp_dir=$(mktemp -d)
        local new_datasets=()
        
        for prefix in "${datasets[@]}"; do
            local base=$(basename "$prefix")
            plink --bfile "$prefix" \
                --exclude "${output_prefix}-merge.missnp" \
                --make-bed \
                --out "${temp_dir}/${base}" \
                >> "$LOG_FILE" 2>&1
            new_datasets+=("${temp_dir}/${base}")
        done
        
        # Retry merge
        rm -f "$merge_list"
        for ((i=1; i<${#new_datasets[@]}; i++)); do
            echo "${new_datasets[$i]}" >> "$merge_list"
        done
        
        plink --bfile "${new_datasets[0]}" \
            --merge-list "$merge_list" \
            --make-bed \
            --out "$output_prefix" \
            --allow-no-sex \
            >> "$LOG_FILE" 2>&1
        
        rm -rf "$temp_dir"
    fi
    
    # Verify merge
    if [[ -f "${output_prefix}.bed" ]]; then
        local n_variants=$(wc -l < "${output_prefix}.bim")
        local n_samples=$(wc -l < "${output_prefix}.fam")
        echo "Merged dataset: ${n_variants} variants, ${n_samples} samples" | tee -a "$LOG_FILE"
        return 0
    else
        echo "ERROR: Merge failed" | tee -a "$LOG_FILE"
        return 1
    fi
}

#-------------------------------------------------------------------------------
# Function: update_fam_with_population
# Description: Add population information to FAM file
#-------------------------------------------------------------------------------
update_fam_with_population() {
    local fam_file="$1"
    local metadata_file="${REF_DIR}/sample_metadata.tsv"
    
    if [[ ! -f "$metadata_file" ]]; then
        echo "Warning: Metadata file not found, skipping population annotation" | tee -a "$LOG_FILE"
        return 0
    fi
    
    echo "Adding population information to FAM file..." | tee -a "$LOG_FILE"
    
    # Create backup
    cp "$fam_file" "${fam_file}.bak"
    
    # Note: PLINK FAM file has limited fields, population info stored separately
    # Create a sample-population mapping file for downstream analysis
    local pop_file="${fam_file%.fam}.pop"
    
    awk -v metadata="$metadata_file" '
    BEGIN {
        while ((getline line < metadata) > 0) {
            split(line, a, "\t");
            pop[a[1]] = a[2];
        }
    }
    {
        sample = $2;
        if (sample in pop) {
            print pop[sample];
        } else {
            print "USER";  # Mark user samples
        }
    }' "$fam_file" > "$pop_file"
    
    echo "Population file created: $pop_file" | tee -a "$LOG_FILE"
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
    
    # Check user data
    echo "=== Checking Input Files ===" | tee -a "$LOG_FILE"
    check_plink_files "$USER_PREFIX" "User data" || exit 1
    
    # Convert reference VCFs to PLINK if needed
    echo "" | tee -a "$LOG_FILE"
    echo "=== Preparing Reference Data ===" | tee -a "$LOG_FILE"
    
    local kg_prefix="${PROCESSED_DIR}/1kg_euro"
    local hgdp_prefix="${PROCESSED_DIR}/hgdp_euro"
    
    if [[ -f "${RAW_DIR}/1kg_euro_merged.vcf.gz" ]]; then
        convert_vcf_to_plink "${RAW_DIR}/1kg_euro_merged.vcf.gz" "$kg_prefix"
    fi
    
    if [[ -f "${RAW_DIR}/hgdp_euro.vcf.gz" ]]; then
        convert_vcf_to_plink "${RAW_DIR}/hgdp_euro.vcf.gz" "$hgdp_prefix"
    fi
    
    # Prepare datasets for merging
    echo "" | tee -a "$LOG_FILE"
    echo "=== Merging Datasets ===" | tee -a "$LOG_FILE"
    
    local datasets_to_merge=()
    datasets_to_merge+=("$USER_PREFIX")
    
    if [[ -f "${kg_prefix}.bed" ]]; then
        handle_strand_flips "$USER_PREFIX" "$kg_prefix" "${PROCESSED_DIR}/1kg_aligned"
        datasets_to_merge+=("${PROCESSED_DIR}/1kg_aligned")
    fi
    
    if [[ -f "${hgdp_prefix}.bed" ]]; then
        handle_strand_flips "$USER_PREFIX" "$hgdp_prefix" "${PROCESSED_DIR}/hgdp_aligned"
        datasets_to_merge+=("${PROCESSED_DIR}/hgdp_aligned")
    fi
    
    # Merge all datasets
    local merged_prefix="${PROCESSED_DIR}/merged_all"
    merge_datasets "$merged_prefix" "${datasets_to_merge[@]}"
    
    # Add population information
    update_fam_with_population "${merged_prefix}.fam"
    
    echo "" | tee -a "$LOG_FILE"
    echo "============================================" | tee -a "$LOG_FILE"
    echo "Data merge completed successfully!" | tee -a "$LOG_FILE"
    echo "Output: ${merged_prefix}.bed/bim/fam" | tee -a "$LOG_FILE"
    echo "Finished: $(date)" | tee -a "$LOG_FILE"
    echo "============================================" | tee -a "$LOG_FILE"
}

# Run main function
main "$@"
