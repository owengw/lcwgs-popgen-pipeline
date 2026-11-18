#!/bin/bash
#SBATCH --job-name=admixture
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=48:00:00
#SBATCH --array=1-200
#SBATCH --output=logs/admixture_%A_%a.log

# ADMIXTURE analysis with K=1-10, 20 reps each
# Array jobs: 10 K values × 20 reps = 200 jobs

source ~/.bash_profile
conda activate /mnt/parscratch/users/bi4og/conda_envs/owenspopgen

set -eu

# Variables
OUT_DIR=${OUT_DIR:-physalia/angsd_results}
ADMIXTURE_DIR="$OUT_DIR/admixture"
BEAGLE_FILE="$OUT_DIR/pop_map_ALL.beagle.gz"
PLINK_PREFIX="$ADMIXTURE_DIR/pop_map_ALL"

mkdir -p "$ADMIXTURE_DIR"

# Calculate K and rep from array task ID
K=$(( (SLURM_ARRAY_TASK_ID - 1) / 20 + 1 ))
REP=$(( (SLURM_ARRAY_TASK_ID - 1) % 20 + 1 ))

echo "=================================================================="
echo "ADMIXTURE Analysis: K=$K, Rep=$REP"
echo "=================================================================="

# Convert beagle to PLINK format (only for first job)
if [[ $SLURM_ARRAY_TASK_ID -eq 1 ]]; then
    echo "Converting beagle to PLINK format..."
    
    if [[ ! -f "$BEAGLE_FILE" ]]; then
        echo "ERROR: Beagle file not found: $BEAGLE_FILE"
        exit 1
    fi
    
    # Convert using PLINK via intermediate PED format
    Rscript --vanilla - "$BEAGLE_FILE" "$PLINK_PREFIX" <<'REOF'
suppressPackageStartupMessages(library(tidyverse))

args <- commandArgs(trailingOnly = TRUE)
beagle_file <- args[1]
out_prefix <- args[2]

cat("Reading beagle file...\n")
beagle <- read_tsv(beagle_file, col_types = cols(), show_col_types = FALSE)

n_sites <- nrow(beagle)
n_ind <- (ncol(beagle) - 3) / 3

cat("Sites:", n_sites, "Individuals:", n_ind, "\n")

# Create PED file data
ped_matrix <- matrix(nrow = n_ind, ncol = 6 + 2 * n_sites)

for (i in 1:n_ind) {
  # PLINK format: FID IID PID MID SEX PHENO + genotypes
  ped_matrix[i, 1] <- paste0("FAM", i)
  ped_matrix[i, 2] <- paste0("IND", i)
  ped_matrix[i, 3:6] <- c("0", "0", "0", "-9")
  
  for (site in 1:n_sites) {
    # Get genotype likelihoods for this individual at this site
    gl_col_start <- 3 + (i - 1) * 3 + 1
    gl_vals <- as.numeric(beagle[site, gl_col_start:(gl_col_start + 2)])
    
    # Get alleles and convert numeric codes to nucleotides
    allele_ref <- as.character(beagle$allele1[site])
    allele_alt <- as.character(beagle$allele2[site])
    
    # Convert numeric allele codes (0,1,2) to nucleotides (A,C,G,T)
    # This avoids PLINK confusion with 0 as missing
    allele_map <- c("0" = "A", "1" = "C", "2" = "G", "3" = "T")
    
    if (allele_ref %in% names(allele_map)) {
      allele_ref <- allele_map[allele_ref]
    }
    if (allele_alt %in% names(allele_map)) {
      allele_alt <- allele_map[allele_alt]
    }
    
    # Check for missing data: all three GLs are equal or any is NA
    is_missing <- (length(unique(gl_vals)) == 1) || any(is.na(gl_vals))
    
    if (is_missing) {
      # Mark as missing in PLINK format (0 0)
      ped_matrix[i, 6 + 2 * site - 1] <- "0"
      ped_matrix[i, 6 + 2 * site] <- "0"
    } else {
      # Call genotype based on highest likelihood
      max_idx <- which.max(gl_vals)
      
      if (max_idx == 1) {
        # Homozygous reference (0/0)
        ped_matrix[i, 6 + 2 * site - 1] <- allele_ref
        ped_matrix[i, 6 + 2 * site] <- allele_ref
      } else if (max_idx == 2) {
        # Heterozygous (0/1)
        ped_matrix[i, 6 + 2 * site - 1] <- allele_ref
        ped_matrix[i, 6 + 2 * site] <- allele_alt
      } else {
        # Homozygous alternate (1/1)
        ped_matrix[i, 6 + 2 * site - 1] <- allele_alt
        ped_matrix[i, 6 + 2 * site] <- allele_alt
      }
    }
  }
  
  if (i %% 50 == 0) cat("  Processed", i, "/", n_ind, "individuals\n")
}

# Write PED file
cat("Writing PED file...\n")
ped_file <- paste0(out_prefix, ".ped")
write.table(ped_matrix, ped_file, 
            quote = FALSE, row.names = FALSE, col.names = FALSE, sep = " ")

# Create MAP file
cat("Writing MAP file...\n")
map_df <- data.frame(
  chr = str_extract(beagle$marker, "^[^_]+"),
  snp = beagle$marker,
  cm = 0,
  bp = as.integer(str_extract(beagle$marker, "[0-9]+$"))
)

map_file <- paste0(out_prefix, ".map")
write.table(map_df, map_file,
            quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "\t")

cat("Conversion complete!\n")
cat("  PED file:", ped_file, "\n")
cat("  MAP file:", map_file, "\n")
REOF

    # Convert to binary PLINK format
    echo "Converting to binary PLINK format..."
    plink --file $PLINK_PREFIX \
          --make-bed \
          --out ${PLINK_PREFIX}-temporary \
          --allow-extra-chr \
          --allow-no-sex \
          --geno 0.1 \
          --mind 0.5 2>&1 | tee ${PLINK_PREFIX}_conversion.log
    
    # Recode chromosome names to integers for ADMIXTURE
    # ADMIXTURE requires integer chromosome codes
    echo "Recoding chromosomes to integers for ADMIXTURE..."
    
    # Get unique chromosome names and assign integer codes
    awk '{print $1}' ${PLINK_PREFIX}-temporary.bim | sort -u > ${ADMIXTURE_DIR}/chr_list.txt
    
    # Create a mapping file: old_chr -> new_chr (1, 2, 3, ...)
    awk '{print $0"\t"NR}' ${ADMIXTURE_DIR}/chr_list.txt > ${ADMIXTURE_DIR}/chr_map.txt
    
    # Recode the BIM file
    awk 'NR==FNR {map[$1]=$2; next} {$1=map[$1]; print}' \
        ${ADMIXTURE_DIR}/chr_map.txt \
        ${PLINK_PREFIX}-temporary.bim > ${PLINK_PREFIX}.bim
    
    # Copy the other files
    cp ${PLINK_PREFIX}-temporary.bed ${PLINK_PREFIX}.bed
    cp ${PLINK_PREFIX}-temporary.fam ${PLINK_PREFIX}.fam
    
    echo "Chromosome recoding complete"
    echo "First 5 SNPs after recoding:"
    head -5 ${PLINK_PREFIX}.bim
    
    # Wait for PLINK to finish writing files (it creates temporary files first)
    echo "Waiting for PLINK to complete file operations..."
    MAX_WAIT=60
    WAITED=0
    while [[ ! -f "${PLINK_PREFIX}.bed" ]] && [[ $WAITED -lt $MAX_WAIT ]]; do
        sleep 1
        WAITED=$((WAITED + 1))
    done
    
    if [[ ! -f "${PLINK_PREFIX}.bed" ]]; then
        echo "ERROR: PLINK conversion failed - .bed file not found after ${MAX_WAIT}s"
        echo "Files in admixture directory:"
        ls -lh ${ADMIXTURE_DIR}/*.bed ${ADMIXTURE_DIR}/*.bim ${ADMIXTURE_DIR}/*.fam 2>/dev/null || echo "No PLINK binary files found"
        echo ""
        echo "Conversion log:"
        tail -20 ${PLINK_PREFIX}_conversion.log
        exit 1
    fi
    
    # Additional sync to ensure filesystem is caught up
    sync
    sleep 2
    
    # Clean up text format and temporary files
    rm -f ${PLINK_PREFIX}.ped ${PLINK_PREFIX}.map
    rm -f ${PLINK_PREFIX}-temporary.*
    
    echo "Conversion complete!"
    
    # Extract QC stats and save to metadata file
    N_SNPS=$(wc -l < ${PLINK_PREFIX}.bim)
    N_IND=$(wc -l < ${PLINK_PREFIX}.fam)
    
    echo "SNPs retained: $N_SNPS"
    echo "Individuals retained: $N_IND"
    
    # Save QC metadata for report
    cat > ${ADMIXTURE_DIR}/conversion_metadata.txt <<METADATA
# ADMIXTURE Input Data QC Summary
# Generated: $(date)

Original SNPs: 5460
SNPs after QC (--geno 0.1): $N_SNPS
SNPs removed: $((5460 - N_SNPS))

Original individuals: 198
Individuals after QC (--mind 0.5): $N_IND
Individuals removed: $((198 - N_IND))

Mean genotyping rate: $(grep "Total genotyping rate" ${PLINK_PREFIX}_conversion.log | awk '{print $NF}')

QC filters applied:
  --geno 0.1   : Remove SNPs with >10% missing data
  --mind 0.5   : Remove individuals with >50% missing data
  --allow-extra-chr : Allow non-standard chromosome names
  Chromosome recoding: Contig names recoded to integers 1-N for ADMIXTURE compatibility

Input file: pop_map_ALL.beagle.gz
LD-pruned SNPs: Yes (from PCA analysis)
METADATA
    
    echo "QC metadata saved to: ${ADMIXTURE_DIR}/conversion_metadata.txt"
    
    # Create a marker file to indicate conversion is complete
    touch ${PLINK_PREFIX}.conversion_complete
fi

# Wait for conversion to complete (for jobs 2-200)
if [[ $SLURM_ARRAY_TASK_ID -gt 1 ]]; then
    echo "Waiting for PLINK conversion to complete..."
    while [[ ! -f "${PLINK_PREFIX}.conversion_complete" ]]; do
        sleep 10
    done
    echo "Conversion complete, proceeding..."
fi

# Create working directory for this K and rep
WORK_DIR="$ADMIXTURE_DIR/K${K}_rep${REP}"
mkdir -p "$WORK_DIR"

# Copy PLINK files to working directory with retry logic
echo "Copying PLINK files to working directory..."
MAX_RETRIES=5
RETRY=0

while [[ $RETRY -lt $MAX_RETRIES ]]; do
    # Use full paths since we're in the admixture directory
    if [[ -f "${PLINK_PREFIX}.bed" ]] && [[ -f "${PLINK_PREFIX}.bim" ]] && [[ -f "${PLINK_PREFIX}.fam" ]]; then
        cp ${PLINK_PREFIX}.bed ${WORK_DIR}/pop_map_ALL.bed && \
        cp ${PLINK_PREFIX}.bim ${WORK_DIR}/pop_map_ALL.bim && \
        cp ${PLINK_PREFIX}.fam ${WORK_DIR}/pop_map_ALL.fam
        
        if [[ $? -eq 0 ]]; then
            echo "Files copied successfully"
            break
        fi
    fi
    
    RETRY=$((RETRY + 1))
    echo "Copy attempt $RETRY failed, waiting 5 seconds..."
    sleep 5
done

if [[ $RETRY -eq $MAX_RETRIES ]]; then
    echo "ERROR: Failed to copy PLINK files after $MAX_RETRIES attempts"
    echo "Source files status:"
    ls -lh ${PLINK_PREFIX}.bed ${PLINK_PREFIX}.bim ${PLINK_PREFIX}.fam 2>&1
    exit 1
fi

# Now change to working directory
cd "$WORK_DIR"

# Run ADMIXTURE with different random seed for each replicate
SEED=$((K * 1000 + REP))

echo "Running ADMIXTURE (K=$K, seed=$SEED)..."
admixture -s $SEED --cv pop_map_ALL.bed $K | tee admixture_K${K}_rep${REP}.log

# Check completion
if [[ -f "pop_map_ALL.${K}.Q" ]]; then
    # Rename output files
    mv pop_map_ALL.${K}.Q admixture_K${K}_rep${REP}.Q
    mv pop_map_ALL.${K}.P admixture_K${K}_rep${REP}.P
    
    # Extract CV error
    CV_ERROR=$(grep "CV error" admixture_K${K}_rep${REP}.log | awk '{print $4}')
    echo "K=$K Rep=$REP CV_error=$CV_ERROR" >> ../cv_errors.txt
    
    echo "ADMIXTURE completed successfully"
    echo "  CV error: $CV_ERROR"
else
    echo "ERROR: ADMIXTURE failed"
    exit 1
fi

# Clean up
rm -f pop_map_ALL.bed pop_map_ALL.bim pop_map_ALL.fam

echo "Complete: K=$K, Rep=$REP"