#!/bin/bash
#SBATCH --job-name=lcmlkin
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=12:00:00
#SBATCH --output=logs/lcmlkin_%j.log

# =============================================================================
# lcMLkin v2.1 relatedness estimation
# Uses VCF from smcpp pipeline (pllong_main_pops.vcf.gz)
# Runs per population: CAI, NS, SB, SI
#
# lcMLkin v2.1 preferred over ngsRelate KING for highly inbred populations —
# uses a likelihood approach that accounts for inbreeding directly rather than
# assuming HWE as KING requires.
#
# Output columns: Z0ag=K0, Z1ag=K1, Z2ag=K2, PI_HATag=r
# Relationship thresholds:
#   Duplicate:    K0~0.0, r~1.0
#   1st degree:   K0~0.0, r~0.5  (parent-offspring)
#   Full sibling: K0~0.25, K1~0.5, K2~0.25, r~0.5
#   2nd degree:   K0~0.5, r~0.25
#   3rd degree:   K0~0.75, r~0.125
#   Unrelated:    K0~1.0, r~0.0
#
# Key input preparation steps:
#   - VCF must have . SNP IDs (lcMLkin constructs CHROM_POS internally)
#   - Plink bim must have CHROM_POS IDs to match lcMLkin's internal IDs
#   - Freq file must be pre-generated so lcMLkin finds it and skips its
#     internal plink freq call (which fails on non-standard chromosomes)
#   - lcMLkin must be run from the output directory so basename paths resolve
# =============================================================================

source ~/.bash_profile
module load Python/3.10.8-GCCcore-12.2.0
module load R/4.4.1-foss-2022b 2>/dev/null || true
export PYTHONPATH=/users/bi4og/.local/lib/python3.10/site-packages:$PYTHONPATH

set -euo pipefail
set -x

# =============================================================================
# Variables
# =============================================================================
OUT_DIR=${OUT_DIR:-physalia/angsd_results}
SMCPP_DIR="physalia/smcpp"
VCF_ALL="$SMCPP_DIR/pllong_main_pops.vcf.gz"
LCMLKIN_DIR="/mnt/parscratch/users/bi4og/lcMLkin-v2.1"
LCMLKIN_SCRIPT="$LCMLKIN_DIR/lcmlkinv2.py"
LCMLKIN_OUT="$OUT_DIR/lcmlkin"
POP_LIST_DIR="$OUT_DIR/pop_map_lists"
PLINK="/mnt/community/Genomics/bin/plink"
BCFTOOLS=$(which bcftools 2>/dev/null || echo "/mnt/parscratch/users/bi4og/conda_envs/owenspopgen/bin/bcftools")
EXCLUDE="o_merged\|T_merged\|SI43_merged\|SI45_merged\|SI83_merged"
R_LIBS=/mnt/parscratch/users/bi4og/R_libs
THREADS=$SLURM_CPUS_PER_TASK

POPULATIONS=(CAI NS SB SI)

mkdir -p "$LCMLKIN_OUT"
mkdir -p logs

# =============================================================================
# Step 1 — Validate inputs
# =============================================================================
echo "=========================================="
echo "Validating inputs"
echo "=========================================="
[[ ! -f "$VCF_ALL" ]]        && echo "ERROR: VCF not found: $VCF_ALL"               && exit 1
[[ ! -f "$PLINK" ]]          && echo "ERROR: plink not found: $PLINK"                && exit 1
[[ ! -f "$LCMLKIN_SCRIPT" ]] && echo "ERROR: lcMLkin not found: $LCMLKIN_SCRIPT"     && exit 1
echo "  VCF, plink, lcMLkin — OK"

python3 -c "import numpy, pandas, scipy; print('  Python packages OK')" || {
    echo "ERROR: Required Python packages not available"
    exit 1
}

# =============================================================================
# Step 2 — Build clean sample ID mapping and rename VCF
# VCF sample names are full BAM paths — extract clean IDs
# VCF must keep . SNP IDs (lcMLkin constructs CHROM_POS internally)
# =============================================================================
echo "=========================================="
echo "Preparing renamed VCF"
echo "=========================================="

SAMPLE_MAP="$LCMLKIN_OUT/all_sample_ids.txt"
RENAMED_VCF="$LCMLKIN_OUT/main_pops_renamed.vcf.gz"

if [[ ! -f "$RENAMED_VCF" ]]; then
    echo "  Building sample name mapping..."
    "$BCFTOOLS" query -l "$VCF_ALL" | \
        sed 's|.*/||' | \
        sed 's/_merged_pe_bt2.*$//' | \
        sed 's/_merged_se_bt2.*$//' > "$SAMPLE_MAP"

    N_SAMPLES=$(wc -l < "$SAMPLE_MAP")
    echo "  Samples mapped: $N_SAMPLES"

    echo "  Renaming VCF samples..."
    "$BCFTOOLS" reheader --samples "$SAMPLE_MAP" "$VCF_ALL" | \
        "$BCFTOOLS" view --output-type z --output "$RENAMED_VCF"
    "$BCFTOOLS" index -t "$RENAMED_VCF"
    echo "  Renamed VCF ready"
else
    echo "  Renamed VCF already exists — skipping"
fi

# =============================================================================
# Step 3 — Per-population: prepare inputs and run lcMLkin
# =============================================================================
for POP in "${POPULATIONS[@]}"; do
    echo "=========================================="
    echo "Processing: $POP"
    echo "=========================================="

    POP_OUT="$LCMLKIN_OUT/$POP"
    mkdir -p "$POP_OUT"

    FULL_LIST="$POP_LIST_DIR/pop_map_${POP}.list"
    if [[ ! -f "$FULL_LIST" ]]; then
        echo "  WARNING: BAM list not found: $FULL_LIST — skipping"
        continue
    fi

    # --- 3a: Build clean sample list ---
    SAMPLE_LIST="$POP_OUT/${POP}_samples.txt"
    grep -v "$EXCLUDE" "$FULL_LIST" | \
        sed 's|.*/||' | \
        sed 's/_merged_pe_bt2.*$//' | \
        sed 's/_merged_se_bt2.*$//' > "$SAMPLE_LIST"

    N_POP=$(wc -l < "$SAMPLE_LIST")
    echo "  Samples (excl. duplicates/failures): $N_POP"
    [[ $N_POP -lt 2 ]] && echo "  WARNING: fewer than 2 samples — skipping" && continue

    # Verify samples exist in VCF
    N_MATCHED=$("$BCFTOOLS" query -l "$RENAMED_VCF" | grep -Fxf "$SAMPLE_LIST" | wc -l)
    echo "  Samples matched in VCF: $N_MATCHED / $N_POP"
    if [[ $N_MATCHED -lt 2 ]]; then
        echo "  ERROR: fewer than 2 samples matched in VCF"
        continue
    fi
    # --- 3b: Use ANGSD VCF directly (already per-population) ---
    POP_VCF="$POP_OUT/${POP}_angsd.vcf.gz"
    if [[ ! -f "$POP_VCF" ]]; then
        echo "  ERROR: ANGSD VCF not found: $POP_VCF"
        echo "  Run p09a_angsd_vcf.sh first"
        continue
    fi
    N_SITES=$("$BCFTOOLS" stats "$POP_VCF" | grep "^SN.*number of SNPs" | awk '{print $NF}')
    N_SAMP=$("$BCFTOOLS" query -l "$POP_VCF" | wc -l)
    echo "  VCF: $N_SAMP samples, $N_SITES SNPs"


    # --- 3b.2: Filter to biallelic SNPs with MAF>0.05 and thin to ~50k ---
    FILTERED_VCF="$POP_OUT/${POP}_filtered.vcf.gz"
    if [[ ! -f "$FILTERED_VCF" ]]; then
        echo "  Filtering VCF (biallelic, MAF>0.05, thin to 50k)..."
        "$BCFTOOLS" view \
            --min-af 0.05:minor \
            --max-alleles 2 \
            --min-alleles 2 \
            "$POP_VCF" 2>/dev/null | \
        "$BCFTOOLS" view \
            --output-type z \
            --output "$FILTERED_VCF"
        "$BCFTOOLS" index -t "$FILTERED_VCF"
        N_FILT=$("$BCFTOOLS" stats "$FILTERED_VCF" 2>/dev/null | \
                 grep "^SN.*number of SNPs" | awk '{print $NF}')
        echo "  Filtered VCF: $N_FILT SNPs"
        # If still >100k SNPs, randomly subsample to 100k
        if [[ $N_FILT -gt 50000 ]]; then
            echo "  Subsampling to 100k SNPs..."
            "$BCFTOOLS" view "$FILTERED_VCF" 2>/dev/null | \
            awk 'BEGIN{srand(42)} /^#/{print} !/^#/ && rand()<0.01{print}' | \
            "$BCFTOOLS" view --output-type z --output "${FILTERED_VCF%.vcf.gz}_sub.vcf.gz"
            "$BCFTOOLS" index -t "${FILTERED_VCF%.vcf.gz}_sub.vcf.gz"
            mv "${FILTERED_VCF%.vcf.gz}_sub.vcf.gz" "$FILTERED_VCF"
            mv "${FILTERED_VCF%.vcf.gz}_sub.vcf.gz.tbi" "${FILTERED_VCF%.vcf.gz}.vcf.gz.tbi" 2>/dev/null || true
            N_FILT=$("$BCFTOOLS" stats "$FILTERED_VCF" 2>/dev/null | \
                     grep "^SN.*number of SNPs" | awk '{print $NF}')
            echo "  After subsampling: $N_FILT SNPs"
        fi
    else
        echo "  Filtered VCF already exists — skipping"
        N_FILT=$("$BCFTOOLS" stats "$FILTERED_VCF" 2>/dev/null | \
                 grep "^SN.*number of SNPs" | awk '{print $NF}')
        echo "  Filtered VCF: $N_FILT SNPs"
    fi
    # Strip PL field — lcMLkin fails on PL=. for missing genotypes; GL is sufficient
    echo "  Removing PL field from VCF..."
    "$BCFTOOLS" annotate \
        --remove FORMAT/PL,FORMAT/GQ \
        --output-type z \
        --output "${FILTERED_VCF%.vcf.gz}_nopl.vcf.gz" \
        "$FILTERED_VCF"
    "$BCFTOOLS" index -t "${FILTERED_VCF%.vcf.gz}_nopl.vcf.gz"
    mv "${FILTERED_VCF%.vcf.gz}_nopl.vcf.gz" "$FILTERED_VCF"
    mv "${FILTERED_VCF%.vcf.gz}_nopl.vcf.gz.tbi" "${FILTERED_VCF%.vcf.gz}.vcf.gz.tbi" 2>/dev/null || \
        "$BCFTOOLS" index -t "$FILTERED_VCF"
    POP_VCF="$FILTERED_VCF"
















    # --- 3c: Convert to plink with CHROM_POS IDs in bim ---
    PLINK_PREFIX="$POP_OUT/${POP}_plink"
    if [[ ! -f "${PLINK_PREFIX}.bed" ]]; then
        echo "  Converting to plink binary format..."
        "$PLINK" \
            --vcf "$POP_VCF" \
            --make-bed \
            --allow-extra-chr \
            --double-id \
            --out "$PLINK_PREFIX" \
            > "$POP_OUT/${POP}_plink.log" 2>&1

        if [[ ! -f "${PLINK_PREFIX}.bed" ]]; then
            echo "  ERROR: plink conversion failed"
            tail -10 "$POP_OUT/${POP}_plink.log"
            continue
        fi

        # Add CHROM_POS IDs to bim so --extract matches lcMLkin's internal IDs
        echo "  Adding CHROM_POS IDs to bim file..."
        awk '{$2=$1"_"$4; print}' OFS='\t' "${PLINK_PREFIX}.bim" > "${PLINK_PREFIX}_ids.bim"
        mv "${PLINK_PREFIX}.bim" "${PLINK_PREFIX}.bim.bak"
        mv "${PLINK_PREFIX}_ids.bim" "${PLINK_PREFIX}.bim"
        echo "  Plink files ready: $(wc -l < ${PLINK_PREFIX}.bim) variants"
    else
        echo "  Plink files already exist — checking bim IDs..."
        # Ensure bim has CHROM_POS IDs (fix if it has . IDs)
        if awk 'NR==1{exit ($2==".")?0:1}' "${PLINK_PREFIX}.bim"; then
            echo "  Bim has . IDs — adding CHROM_POS IDs..."
            awk '{$2=$1"_"$4; print}' OFS='\t' "${PLINK_PREFIX}.bim" > "${PLINK_PREFIX}_ids.bim"
            mv "${PLINK_PREFIX}.bim" "${PLINK_PREFIX}.bim.bak"
            mv "${PLINK_PREFIX}_ids.bim" "${PLINK_PREFIX}.bim"
        else
            echo "  Bim already has CHROM_POS IDs — OK"
        fi
    fi

    # --- 3d: Pre-generate freq file ---
    # lcMLkin checks for this file first; pre-generating avoids its internal
    # plink freq call which uses --extract with CHROM_POS IDs that would
    # fail against a . ID bim
    echo "  Pre-generating allele frequency file..."
    "$PLINK" \
        --bfile "$PLINK_PREFIX" \
        --freq \
        --allow-extra-chr \
        --out "$PLINK_PREFIX" \
        > "$POP_OUT/${POP}_plink_freq.log" 2>&1

    if [[ ! -f "${PLINK_PREFIX}.frq" ]]; then
        echo "  ERROR: freq file not generated"
        tail -5 "$POP_OUT/${POP}_plink_freq.log"
        continue
    fi
    echo "  Freq file ready: $(wc -l < ${PLINK_PREFIX}.frq) variants"

    # --- 3e: Run lcMLkin from output directory ---
    # Must cd to POP_OUT so os.path.basename() paths resolve correctly
    RESULT_FILE="$POP_OUT/${POP}_lcmlkin.txt"
    if [[ -f "$RESULT_FILE" ]]; then
        N_PAIRS=$(wc -l < "$RESULT_FILE")
        echo "  lcMLkin results already exist ($N_PAIRS pairs) — skipping"
        continue
    fi

    echo "  Running lcMLkin from: $POP_OUT"
    cd "$POP_OUT"
    python3 "$LCMLKIN_SCRIPT" \
        -v "${POP}_filtered.vcf.gz" \
        -p "${POP}_plink" \
        -t "$THREADS" \
        -o "${POP}_lcmlkin" \
        > "${POP}_lcmlkin.log" 2>&1

    cd /mnt/parscratch/users/bi4og

    if [[ ! -f "$RESULT_FILE" ]]; then
        echo "  ERROR: No result file produced — check $POP_OUT/${POP}_lcmlkin.log"
        tail -10 "$POP_OUT/${POP}_lcmlkin.log"
        continue
    fi

    N_PAIRS=$(wc -l < "$RESULT_FILE")
    echo "  Done: $POP — $N_PAIRS pairs"
done

# =============================================================================
# Step 4 — Collate results across populations
# =============================================================================
echo ""
echo "=========================================="
echo "Collating results"
echo "=========================================="

export R_LIBS_USER="$R_LIBS"

cat > /tmp/collate_lcmlkin.R << 'REOF'
pop_order   <- c("CAI", "NS", "SB", "SI")
lcmlkin_dir <- "/mnt/parscratch/users/bi4og/physalia/angsd_results/lcmlkin"

classify_relationship <- function(K0, r) {
    ifelse(r >= 0.45,                              "Duplicate",
    ifelse(K0 <= 0.1  & r >= 0.4,                 "1st degree",
    ifelse(K0 >= 0.15 & K0 <= 0.35 & r >= 0.35,   "Full sibling",
    ifelse(K0 >= 0.4  & K0 <= 0.6,                "2nd degree",
    ifelse(K0 >= 0.65 & K0 <= 0.85,               "3rd degree",
                                                   "Unrelated")))))
}

all_results <- list()
for (pop in pop_order) {
    result_file <- file.path(lcmlkin_dir, pop,
                             paste0(pop, "_lcmlkin.txt"))
    if (!file.exists(result_file)) {
        cat("WARNING: No result file for", pop, "\n")
        next
    }
    res <- tryCatch(
        read.table(result_file, header = TRUE, stringsAsFactors = FALSE),
        error = function(e) { cat("ERROR:", conditionMessage(e), "\n"); NULL }
    )
    if (is.null(res)) next
    cat(pop, "— rows:", nrow(res),
        "cols:", paste(names(res), collapse=", "), "\n")
    res$population <- pop
    all_results[[pop]] <- res
}

if (length(all_results) == 0) {
    cat("No results to collate\n"); quit(status=1)
}

combined <- do.call(rbind, all_results)

if (all(c("Z0ag", "PI_HATag") %in% names(combined))) {
    combined$relationship <- classify_relationship(combined$Z0ag,
                                                   combined$PI_HATag)
    cat("\n=== Relationship class summary ===\n")
    print(table(combined$population, combined$relationship))

    cat("\n=== Close pairs (r >= 0.125) ===\n")
    close <- combined[!is.na(combined$PI_HATag) &
                      combined$PI_HATag >= 0.125, ]
    close <- close[order(close$population, -close$PI_HATag), ]
    print(close[, c("population", "Ind1", "Ind2",
                    "Z0ag", "Z1ag", "Z2ag", "PI_HATag", "relationship")])
}

out_file <- file.path(lcmlkin_dir, "lcmlkin_all_populations.tsv")
write.table(combined, out_file, sep="\t", row.names=FALSE, quote=FALSE)
cat("\nWritten to:", out_file, "\n")
REOF

Rscript --vanilla /tmp/collate_lcmlkin.R

echo ""
echo "=========================================="
echo "lcMLkin complete. Results in: $LCMLKIN_OUT"
echo "=========================================="
