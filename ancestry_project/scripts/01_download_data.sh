#!/bin/bash
#===============================================================================
# Script: 01_download_data.sh
# Description: Download reference data from 1000 Genomes Phase 3 and HGDP
#              for European ancestry analysis
#
# Usage: ./01_download_data.sh
#
# Requirements:
#   - wget
#   - bcftools
#   - ~75GB disk space
#
# Output:
#   - data/raw/1kg_euro.vcf.gz (1000 Genomes European samples)
#   - data/raw/hgdp_euro.vcf.gz (HGDP European samples)
#   - reference/sample_metadata.tsv (population metadata)
#===============================================================================

set -euo pipefail

# Configuration
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_DIR="$(dirname "$SCRIPT_DIR")"
RAW_DIR="${PROJECT_DIR}/data/raw"
REF_DIR="${PROJECT_DIR}/reference"
LOG_FILE="${PROJECT_DIR}/results/01_download.log"

# 1000 Genomes Phase 3 URLs
KG_BASE="ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502"

# HGDP URLs (via NCBI or alternative sources)
HGDP_BASE="https://ngs.sanger.ac.uk/production/hgdp"

# European populations from 1000 Genomes
KG_EUR_POPS=("CEU" "GBR" "FIN" "IBS" "TSI")

# European populations from HGDP
HGDP_EUR_POPS=("French" "Italian" "Russian" "Basque" "Sardinian" "Orcadian")

# Create directories
mkdir -p "$RAW_DIR" "$REF_DIR" "$(dirname "$LOG_FILE")"

# Initialize log
echo "============================================" | tee "$LOG_FILE"
echo "Ancestry Analysis - Data Download" | tee -a "$LOG_FILE"
echo "Started: $(date)" | tee -a "$LOG_FILE"
echo "============================================" | tee -a "$LOG_FILE"

#-------------------------------------------------------------------------------
# Function: download_with_retry
# Description: Download a file with retry logic
#-------------------------------------------------------------------------------
download_with_retry() {
    local url="$1"
    local output="$2"
    local max_retries=3
    local retry=0
    
    while [[ $retry -lt $max_retries ]]; do
        if wget -c -q --show-progress -O "$output" "$url"; then
            echo "Successfully downloaded: $output" | tee -a "$LOG_FILE"
            return 0
        fi
        retry=$((retry + 1))
        echo "Retry $retry/$max_retries for $url" | tee -a "$LOG_FILE"
        sleep 10
    done
    
    echo "ERROR: Failed to download $url after $max_retries attempts" | tee -a "$LOG_FILE"
    return 1
}

#-------------------------------------------------------------------------------
# Download 1000 Genomes Phase 3 data
#-------------------------------------------------------------------------------
download_1kg() {
    echo "" | tee -a "$LOG_FILE"
    echo "=== Downloading 1000 Genomes Phase 3 Data ===" | tee -a "$LOG_FILE"
    
    # Download sample panel (population info)
    local panel_url="${KG_BASE}/integrated_call_samples_v3.20130502.ALL.panel"
    local panel_file="${RAW_DIR}/1kg_panel.txt"
    
    if [[ ! -f "$panel_file" ]]; then
        echo "Downloading sample panel..." | tee -a "$LOG_FILE"
        download_with_retry "$panel_url" "$panel_file"
    else
        echo "Sample panel already exists, skipping." | tee -a "$LOG_FILE"
    fi
    
    # Extract European sample IDs
    local euro_samples="${RAW_DIR}/1kg_euro_samples.txt"
    echo "Extracting European sample IDs..." | tee -a "$LOG_FILE"
    grep -E "CEU|GBR|FIN|IBS|TSI" "$panel_file" | cut -f1 > "$euro_samples"
    local n_samples=$(wc -l < "$euro_samples")
    echo "Found $n_samples European samples" | tee -a "$LOG_FILE"
    
    # Download VCF files for each chromosome
    local merged_vcf="${RAW_DIR}/1kg_euro_merged.vcf.gz"
    
    if [[ ! -f "$merged_vcf" ]]; then
        local vcf_files=()
        
        for chr in {1..22}; do
            local vcf_url="${KG_BASE}/ALL.chr${chr}.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz"
            local vcf_file="${RAW_DIR}/1kg_chr${chr}.vcf.gz"
            local vcf_euro="${RAW_DIR}/1kg_chr${chr}_euro.vcf.gz"
            
            if [[ ! -f "$vcf_euro" ]]; then
                # Download full VCF
                if [[ ! -f "$vcf_file" ]]; then
                    echo "Downloading chromosome $chr..." | tee -a "$LOG_FILE"
                    download_with_retry "$vcf_url" "$vcf_file"
                    download_with_retry "${vcf_url}.tbi" "${vcf_file}.tbi"
                fi
                
                # Subset to European samples
                echo "Extracting European samples from chr$chr..." | tee -a "$LOG_FILE"
                bcftools view -S "$euro_samples" -Oz -o "$vcf_euro" "$vcf_file"
                bcftools index -t "$vcf_euro"
            fi
            
            vcf_files+=("$vcf_euro")
        done
        
        # Merge all chromosomes
        echo "Merging all chromosomes..." | tee -a "$LOG_FILE"
        bcftools concat -Oz -o "$merged_vcf" "${vcf_files[@]}"
        bcftools index -t "$merged_vcf"
    else
        echo "Merged 1KG European VCF already exists, skipping." | tee -a "$LOG_FILE"
    fi
    
    echo "1000 Genomes download complete!" | tee -a "$LOG_FILE"
}

#-------------------------------------------------------------------------------
# Download HGDP data
#-------------------------------------------------------------------------------
download_hgdp() {
    echo "" | tee -a "$LOG_FILE"
    echo "=== Downloading HGDP Data ===" | tee -a "$LOG_FILE"
    
    local hgdp_vcf="${RAW_DIR}/hgdp.vcf.gz"
    local hgdp_metadata="${RAW_DIR}/hgdp_metadata.txt"
    local hgdp_euro_vcf="${RAW_DIR}/hgdp_euro.vcf.gz"
    
    # Download HGDP VCF (alternative: use gnomAD HGDP subset)
    if [[ ! -f "$hgdp_vcf" ]]; then
        echo "Downloading HGDP genotypes..." | tee -a "$LOG_FILE"
        # Note: URL may need to be updated based on current data availability
        local hgdp_url="${HGDP_BASE}/hgdp_wgs.20190516.full.chr1-22.vcf.gz"
        download_with_retry "$hgdp_url" "$hgdp_vcf"
        download_with_retry "${hgdp_url}.tbi" "${hgdp_vcf}.tbi"
    fi
    
    # Download metadata
    if [[ ! -f "$hgdp_metadata" ]]; then
        echo "Downloading HGDP metadata..." | tee -a "$LOG_FILE"
        local meta_url="${HGDP_BASE}/hgdp_wgs.20190516.metadata.txt"
        download_with_retry "$meta_url" "$hgdp_metadata"
    fi
    
    # Extract European samples
    if [[ ! -f "$hgdp_euro_vcf" ]]; then
        local euro_samples="${RAW_DIR}/hgdp_euro_samples.txt"
        echo "Extracting European samples from HGDP..." | tee -a "$LOG_FILE"
        
        # Filter for European populations
        grep -E "French|Italian|Russian|Basque|Sardinian|Orcadian" "$hgdp_metadata" | \
            cut -f1 > "$euro_samples"
        
        local n_samples=$(wc -l < "$euro_samples")
        echo "Found $n_samples HGDP European samples" | tee -a "$LOG_FILE"
        
        # Subset VCF
        bcftools view -S "$euro_samples" -Oz -o "$hgdp_euro_vcf" "$hgdp_vcf"
        bcftools index -t "$hgdp_euro_vcf"
    else
        echo "HGDP European VCF already exists, skipping." | tee -a "$LOG_FILE"
    fi
    
    echo "HGDP download complete!" | tee -a "$LOG_FILE"
}

#-------------------------------------------------------------------------------
# Create unified metadata file
#-------------------------------------------------------------------------------
create_metadata() {
    echo "" | tee -a "$LOG_FILE"
    echo "=== Creating Unified Metadata ===" | tee -a "$LOG_FILE"
    
    local metadata_file="${REF_DIR}/sample_metadata.tsv"
    
    # Create header
    echo -e "sample_id\tpopulation\tcountry\tsource" > "$metadata_file"
    
    # Add 1000 Genomes samples
    local kg_panel="${RAW_DIR}/1kg_panel.txt"
    if [[ -f "$kg_panel" ]]; then
        grep -E "CEU|GBR|FIN|IBS|TSI" "$kg_panel" | while read -r line; do
            sample=$(echo "$line" | cut -f1)
            pop=$(echo "$line" | cut -f2)
            
            # Map population to country
            case $pop in
                CEU) country="United States (Utah/Northern European)" ;;
                GBR) country="United Kingdom" ;;
                FIN) country="Finland" ;;
                IBS) country="Spain" ;;
                TSI) country="Italy (Tuscan)" ;;
                *) country="Unknown" ;;
            esac
            
            echo -e "${sample}\t${pop}\t${country}\t1000Genomes"
        done >> "$metadata_file"
    fi
    
    # Add HGDP samples
    local hgdp_meta="${RAW_DIR}/hgdp_metadata.txt"
    if [[ -f "$hgdp_meta" ]]; then
        grep -E "French|Italian|Russian|Basque|Sardinian|Orcadian" "$hgdp_meta" | \
        while read -r line; do
            sample=$(echo "$line" | cut -f1)
            pop=$(echo "$line" | cut -f2)
            
            # Map population to country
            case $pop in
                French) country="France" ;;
                Italian) country="Italy" ;;
                Russian) country="Russia" ;;
                Basque) country="Spain (Basque)" ;;
                Sardinian) country="Italy (Sardinia)" ;;
                Orcadian) country="United Kingdom (Orkney)" ;;
                *) country="Unknown" ;;
            esac
            
            echo -e "${sample}\t${pop}\t${country}\tHGDP"
        done >> "$metadata_file"
    fi
    
    local n_samples=$(tail -n +2 "$metadata_file" | wc -l)
    echo "Created metadata for $n_samples reference samples" | tee -a "$LOG_FILE"
    
    # Create population to country mapping
    local mapping_file="${REF_DIR}/population_country_mapping.tsv"
    cat > "$mapping_file" << 'EOF'
population	country	region
CEU	United States	Northern European ancestry
GBR	United Kingdom	British
FIN	Finland	Finnish
IBS	Spain	Iberian
TSI	Italy	Tuscan
French	France	French
Italian	Italy	North Italian
Russian	Russia	Russian
Basque	Spain	Basque
Sardinian	Italy	Sardinian
Orcadian	United Kingdom	Orcadian
EOF
    
    echo "Created population-country mapping" | tee -a "$LOG_FILE"
}

#-------------------------------------------------------------------------------
# Main execution
#-------------------------------------------------------------------------------
main() {
    echo "" | tee -a "$LOG_FILE"
    echo "Starting reference data download..." | tee -a "$LOG_FILE"
    echo "This may take several hours depending on connection speed." | tee -a "$LOG_FILE"
    echo "Estimated disk space required: ~75GB" | tee -a "$LOG_FILE"
    echo "" | tee -a "$LOG_FILE"
    
    # Check for required tools
    for tool in wget bcftools; do
        if ! command -v "$tool" &> /dev/null; then
            echo "ERROR: Required tool '$tool' not found. Please install it first." | tee -a "$LOG_FILE"
            exit 1
        fi
    done
    
    # Download reference datasets
    download_1kg
    download_hgdp
    
    # Create unified metadata
    create_metadata
    
    echo "" | tee -a "$LOG_FILE"
    echo "============================================" | tee -a "$LOG_FILE"
    echo "Data download completed successfully!" | tee -a "$LOG_FILE"
    echo "Finished: $(date)" | tee -a "$LOG_FILE"
    echo "============================================" | tee -a "$LOG_FILE"
}

# Run main function
main "$@"
