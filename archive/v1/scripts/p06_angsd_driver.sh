#!/bin/bash
#SBATCH --job-name=angsd_driver
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=96:00:00
#SBATCH --output=logs/angsd_driver_%A.log

# Load environment
source ~/.bash_profile
conda activate /mnt/parscratch/users/bi4og/conda_envs/owenspopgen

set -euo pipefail
set -x

# -------------------------------
# User-defined arguments
# -------------------------------
helpFunction() {
    echo ""
    echo "ANGSD Population Genomics Pipeline with Auto-metadata"
    echo "======================================================"
    echo ""
    echo "Usage: $0 -o <out_dir> -r <reference> -b <bam_list> [-a ancestral] [-t threads] [-p min_mapq] [-q min_baseq] [-i min_ind_ratio] [-R resume_step] [-P run_pca]"
    echo ""
    echo "Required arguments:"
    echo "  -o    Output directory"
    echo "  -r    Reference genome FASTA file (must be indexed)"
    echo "  -b    File containing all BAM paths (one per line)"
    echo ""
    echo "Optional arguments:"
    echo "  -a    Ancestral genome FASTA file (default: use reference)"
    echo "  -t    Number of threads (default: 8)"
    echo "  -p    Minimum mapping quality (default: 10)"
    echo "  -q    Minimum base quality (default: 13)"
    echo "  -f    Fold SFS: 1=folded, 0=unfolded (default: 1)"
    echo "  -w    Window size in bp (default: 50000)"
    echo "  -s    Window step in bp (default: 10000)"
    echo "  -i    Minimum individual ratio (default: 0.20)"
    echo "  -l    Number of samples for LD pruning (default: 50)"
    echo "  -R    Resume from step (0=full, 1=skip ANGSD, 2=skip LD, 3=skip SAF, 4=skip theta, 5=skip FST)"
    echo "  -P    Run PCA analysis (0=no, 1=yes, default: 0)"
    echo "  -h    Show this help"
    echo ""
    echo "Pipeline steps:"
    echo "  0. Subsampled ANGSD for LD pruning"
    echo "  1. LD pruning with ngsLD"
    echo "  2. ANGSD SAF generation (essential for theta/FST)"
    echo "  2b. ANGSD PCA analysis (optional, run with -P 1)"
    echo "  3. Theta calculation"
    echo "  4. FST calculation"
    echo "  5. Summary statistics"
    echo ""
    exit 1
}

# -------------------------------
# Defaults
# -------------------------------
OUT_DIR=""
REFERENCE=""
BAM_LIST=""
ANCESTRAL=""
THREADS=8
MIN_MAPQ=10
MIN_Q=13
FOLD_SFS=1
WINDOW_SIZE=50000
WINDOW_STEP=10000
MIN_IND_RATIO=0.20
LD_SAMPLE_SIZE=50
RESUME_FROM=0
RUN_PCA=0  # Default: skip PCA

# -------------------------------
# Parse arguments
# -------------------------------
while getopts "o:r:b:a:t:p:q:f:w:s:i:l:R:P:h" opt; do
    case "$opt" in
        o) OUT_DIR="$OPTARG" ;;
        r) REFERENCE="$OPTARG" ;;
        b) BAM_LIST="$OPTARG" ;;
        a) ANCESTRAL="$OPTARG" ;;
        t) THREADS="$OPTARG" ;;
        p) MIN_MAPQ="$OPTARG" ;;
        q) MIN_Q="$OPTARG" ;;
        f) FOLD_SFS="$OPTARG" ;;
        w) WINDOW_SIZE="$OPTARG" ;;
        s) WINDOW_STEP="$OPTARG" ;;
        i) MIN_IND_RATIO="$OPTARG" ;;
        l) LD_SAMPLE_SIZE="$OPTARG" ;;
        R) RESUME_FROM="$OPTARG" ;;
        P) RUN_PCA="$OPTARG" ;;
        h) helpFunction ;;
        *) helpFunction ;;
    esac
done

# -------------------------------
# Validate required arguments
# -------------------------------
if [[ -z "$OUT_DIR" ]] || [[ -z "$REFERENCE" ]] || [[ -z "$BAM_LIST" ]]; then
    echo "ERROR: Missing required arguments"
    helpFunction
fi

echo "==================================================================="
echo "ANGSD Population Genomics Pipeline with LD Pruning"
echo "==================================================================="
echo "Output directory: $OUT_DIR"
echo "Reference genome: $REFERENCE"
echo "BAM list: $BAM_LIST"
echo "Ancestral genome: ${ANCESTRAL:-None (using reference)}"
echo "Threads: $THREADS"
echo "Min MAPQ: $MIN_MAPQ"
echo "Min base quality: $MIN_Q"
echo "Window size: $WINDOW_SIZE"
echo "Min individual ratio: $MIN_IND_RATIO"
echo "LD sample size: $LD_SAMPLE_SIZE"
echo "Resume from step: $RESUME_FROM"
echo "Run PCA analysis: $RUN_PCA"
echo "==================================================================="

mkdir -p "$OUT_DIR/logs"

# Validate reference and BAMs
[[ ! -f "$REFERENCE" ]] && echo "ERROR: Reference not found" && exit 1
[[ ! -f "$REFERENCE.fai" ]] && echo "ERROR: Reference index not found. Run samtools faidx $REFERENCE" && exit 1
[[ ! -f "$BAM_LIST" ]] && echo "ERROR: BAM list not found" && exit 1
[[ ! -s "$BAM_LIST" ]] && echo "ERROR: BAM list empty" && exit 1

TOTAL_BAMS=$(wc -l < "$BAM_LIST")
echo "Total BAMs in list: $TOTAL_BAMS"

# -------------------------------
# Generate metadata and sample lists (if needed)
# -------------------------------
if [[ $RESUME_FROM -eq 0 ]]; then
    echo ""
    echo "==================================================================="
    echo "Auto-generating metadata from BAM filenames"
    echo "==================================================================="

    METADATA_FILE="$OUT_DIR/metadata.tsv"

    Rscript --vanilla - "$BAM_LIST" "$OUT_DIR" <<'REOF'
library(tidyverse)
args <- commandArgs(trailingOnly = TRUE)
bam_list_file <- args[1]
out_dir <- args[2]
metadata_file <- file.path(out_dir, "metadata.tsv")

bam_paths <- read_lines(bam_list_file)
metadata <- tibble(bam_path = bam_paths) %>%
  mutate(
    filename = basename(bam_path),
    prefix = str_extract(filename, "^[^_]+"),
    number = as.integer(str_extract(prefix, "^\\d+")),
    pop_part = str_extract(prefix, "(?<=-)\\w+"),
    population = str_extract(pop_part, "^[A-Za-z]+") %>% str_remove("o$"),
    site = population,
    year = if_else(number >= 100, 2023, 2024),
    sample_id = prefix
  ) %>%
  select(sample_id, bam_path, population, site, year) %>%
  arrange(population, year, sample_id)

cat("\nMetadata generation summary:\n")
cat("Total samples:", nrow(metadata), "\n")
write_tsv(metadata, metadata_file)
cat("\nMetadata written to:", metadata_file, "\n")
REOF

    echo "? Metadata generated"

    # Generate sample grouping lists
    Rscript --vanilla - "$OUT_DIR" <<'REOF'
library(tidyverse)
args <- commandArgs(trailingOnly = TRUE)
out_dir <- args[1]
metadata_file <- file.path(out_dir, "metadata.tsv")
meta <- read_tsv(metadata_file, col_types = cols(), show_col_types = FALSE)

dir.create(file.path(out_dir, "pop_map_lists"), showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(out_dir, "site_year_map_lists"), showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(out_dir, "year_map_lists"), showWarnings = FALSE, recursive = TRUE)

# Pop maps
pop_groups <- meta %>% group_by(population) %>% summarise(bams = list(bam_path), n = n(), .groups = "drop")
for (i in seq_len(nrow(pop_groups))) {
    writeLines(pop_groups$bams[[i]], file.path(out_dir, "pop_map_lists", paste0("pop_map_", pop_groups$population[i], ".list")))
}
writeLines(meta$bam_path, file.path(out_dir, "pop_map_lists", "pop_map_ALL.list"))

# Site-year maps
site_year_groups <- meta %>% mutate(site_year = paste0(site, year)) %>% group_by(site_year) %>% summarise(bams = list(bam_path), .groups = "drop")
for (i in seq_len(nrow(site_year_groups))) {
    writeLines(site_year_groups$bams[[i]], file.path(out_dir, "site_year_map_lists", paste0(site_year_groups$site_year[i], ".list")))
}

# Year maps
year_groups <- meta %>% group_by(year) %>% summarise(bams = list(bam_path), .groups = "drop")
for (i in seq_len(nrow(year_groups))) {
    writeLines(year_groups$bams[[i]], file.path(out_dir, "year_map_lists", paste0(year_groups$year[i], ".list")))
}
cat("\n? All sample lists generated\n")
REOF
fi

# -------------------------------
# Step 0: Subsampled ANGSD for LD pruning
# -------------------------------
if [[ $RESUME_FROM -le 0 ]]; then
    echo ""
    echo "==================================================================="
    echo "STEP 0: Running ANGSD on subsampled data (for LD pruning)"
    echo "==================================================================="
    
    ALL_SAMPLES_LIST="$OUT_DIR/all_samples.list"
    LD_SUBSAMPLE_LIST="$OUT_DIR/ld_subsample.list"

    cp "$OUT_DIR/pop_map_lists/pop_map_ALL.list" "$ALL_SAMPLES_LIST"
    TOTAL_SAMPLES=$(wc -l < "$ALL_SAMPLES_LIST")
    echo "Total unique samples: $TOTAL_SAMPLES"

    # Stratified LD subsample
    Rscript --vanilla - "$OUT_DIR" "$LD_SAMPLE_SIZE" <<'REOF'
library(tidyverse)
args <- commandArgs(trailingOnly=TRUE)
out_dir <- args[1]; target_n <- as.integer(args[2])
metadata <- read_tsv(file.path(out_dir,"metadata.tsv"), col_types=cols(), show_col_types=FALSE)
set.seed(42)
subsample <- metadata %>% group_by(population) %>% slice_sample(prop=min(1,target_n/nrow(metadata))) %>% ungroup() %>% slice_sample(n=min(target_n,nrow(.)))
writeLines(subsample$bam_path,file.path(out_dir,"ld_subsample.list"))
cat("LD subsample:",nrow(subsample),"samples from",n_distinct(subsample$population),"populations\n")
REOF

    LD_SAMPLE_COUNT=$(wc -l < "$LD_SUBSAMPLE_LIST")
    echo "Using $LD_SAMPLE_COUNT samples for LD detection"

    ALL_ANGSD_JOB=$(sbatch --parsable \
        --job-name="angsd_ld_sub" \
        --cpus-per-task=$THREADS \
        --mem=64G \
        --time=24:00:00 \
        --export=ALL,BAM_LIST="$LD_SUBSAMPLE_LIST",REFERENCE="$REFERENCE",OUT_DIR="$OUT_DIR",ANCESTRAL="$ANCESTRAL",THREADS="$THREADS",MIN_MAPQ="$MIN_MAPQ",MIN_Q="$MIN_Q",FOLD_SFS="$FOLD_SFS",MIN_IND_RATIO="0.5" \
        scripts/p06a_all_samples.sh)
    echo "  Submitted: Job $ALL_ANGSD_JOB"
else
    echo "Skipping Step 0 (subsampled ANGSD)"
    ALL_ANGSD_JOB=""
fi

# -------------------------------
# Step 1: LD pruning
# -------------------------------
if [[ $RESUME_FROM -le 1 ]]; then
    echo ""
    echo "==================================================================="
    echo "STEP 1: LD Pruning with ngsLD"
    echo "==================================================================="

    if [[ -n "$ALL_ANGSD_JOB" ]]; then
        LD_DEP="--dependency=afterok:$ALL_ANGSD_JOB"
    else
        LD_DEP=""
    fi

    LD_PRUNE_JOB=$(sbatch --parsable \
        $LD_DEP \
        --job-name="ld_prune" \
        --cpus-per-task=16 \
        --mem=128G \
        --time=24:00:00 \
        --export=ALL,OUT_DIR="$OUT_DIR",THREADS=16,REFERENCE="$REFERENCE" \
        scripts/p06b_ld.sh)
    echo "  Submitted: Job $LD_PRUNE_JOB"
else
    echo "Skipping Step 1 (LD pruning)"
    LD_PRUNE_JOB=""
fi

# -------------------------------
# Step 2: ANGSD SAF per grouping 
# -------------------------------
if [[ $RESUME_FROM -le 2 ]]; then
    echo ""
    echo "==================================================================="
    echo "STEP 2: Running ANGSD SAF per group (random subset of sites)"
    echo "==================================================================="
    echo "NOTE: SAF generation only - essential for theta/FST"

    declare -A ANGSD_JOBS

    for LIST_DIR in pop_map_lists site_year_map_lists year_map_lists; do
        [[ ! -d "$OUT_DIR/$LIST_DIR" ]] && continue

        echo ""
        echo "Processing $LIST_DIR:"

        for LIST_FILE in "$OUT_DIR/$LIST_DIR"/*.list; do
            [[ ! -f "$LIST_FILE" ]] && continue

            GROUP_NAME=$(basename "$LIST_FILE" .list)
            N_SAMPLES=$(wc -l < "$LIST_FILE")

            # Memory scaling by sample size
            if [[ $N_SAMPLES -ge 150 ]]; then
                CPUS=8; MEM=256G; TIME=96:00:00
            elif [[ $N_SAMPLES -ge 100 ]]; then
                CPUS=8; MEM=256G; TIME=72:00:00
            elif [[ $N_SAMPLES -ge 50 ]]; then
                CPUS=6; MEM=128G; TIME=48:00:00
            else
                CPUS=$THREADS; MEM=64G; TIME=24:00:00
            fi

            # Dependency logic
            if [[ -n "${LD_PRUNE_JOB:-}" ]]; then
                DEP="--dependency=afterok:$LD_PRUNE_JOB"
            else
                DEP=""
            fi

            JOB_ID=$(sbatch --parsable \
                $DEP \
                --job-name="saf_${GROUP_NAME}" \
                --cpus-per-task=$CPUS \
                --mem=$MEM \
                --time=$TIME \
                --export=ALL,BAM_LIST="$LIST_FILE",REFERENCE="$REFERENCE",OUT_DIR="$OUT_DIR",ANCESTRAL="$ANCESTRAL",THREADS="$CPUS",MIN_MAPQ="$MIN_MAPQ",MIN_Q="$MIN_Q",FOLD_SFS="$FOLD_SFS",MIN_IND_RATIO="$MIN_IND_RATIO",USE_RANDOM_SITES=1 \
                scripts/p06c_angsd_subset.sh)

            echo "  $GROUP_NAME ($N_SAMPLES samples): Job $JOB_ID | Mem=$MEM"
            ANGSD_JOBS["$GROUP_NAME"]=$JOB_ID
        done
    done

    NUM_GROUPS=${#ANGSD_JOBS[@]}
    echo ""
    echo "Total SAF jobs submitted: $NUM_GROUPS (includes pop_map_ALL)"
else
    echo "Skipping Step 2 (SAF generation)"
    declare -A ANGSD_JOBS
fi

# -------------------------------
# Step 2b: ANGSD PCA (optional, only for pop_map_ALL)
# -------------------------------
if [[ $RUN_PCA -eq 1 ]]; then
    echo ""
    echo "==================================================================="
    echo "STEP 2b: Running ANGSD PCA (pop_map_ALL only)"
    echo "==================================================================="
    echo "NOTE: PCA analysis with -doMajorMinor 3 for memory efficiency"

    ALL_LIST="$OUT_DIR/pop_map_lists/pop_map_ALL.list"
    N_ALL=$(wc -l < "$ALL_LIST")

    # Dependency logic - can run in parallel with SAF or after LD pruning
    if [[ -n "${LD_PRUNE_JOB:-}" ]]; then
        PCA_DEP="--dependency=afterok:$LD_PRUNE_JOB"
    else
        PCA_DEP=""
    fi

    PCA_JOB=$(sbatch --parsable \
        $PCA_DEP \
        --job-name="pca_ALL" \
        --cpus-per-task=40 \
        --mem=160G \
        --time=96:00:00 \
        --export=ALL,BAM_LIST="$ALL_LIST",REFERENCE="$REFERENCE",OUT_DIR="$OUT_DIR",ANCESTRAL="$ANCESTRAL",THREADS=40,MIN_MAPQ="$MIN_MAPQ",MIN_Q="$MIN_Q",MIN_IND_RATIO="$MIN_IND_RATIO",USE_LD_PRUNED=1 \
        scripts/p06c2_pca.sh)

    echo "  Submitted: PCA Job $PCA_JOB (pop_map_ALL, $N_ALL samples)"
else
    echo ""
    echo "Skipping Step 2b (PCA analysis) - use -P 1 to enable"
fi

# If resuming past step 2, populate group list from existing files
if [[ $RESUME_FROM -gt 2 && ${#ANGSD_JOBS[@]} -eq 0 ]]; then
    echo "Detecting existing SAF files for resume..."
    for LIST_DIR in pop_map_lists site_year_map_lists year_map_lists; do
        [[ ! -d "$OUT_DIR/$LIST_DIR" ]] && continue
        for LIST_FILE in "$OUT_DIR/$LIST_DIR"/*.list; do
            [[ ! -f "$LIST_FILE" ]] && continue
            GROUP_NAME=$(basename "$LIST_FILE" .list)
            # Check if SAF exists
            if [[ -f "$OUT_DIR/${GROUP_NAME}.saf.idx" ]]; then
                ANGSD_JOBS["$GROUP_NAME"]="completed"
                echo "  Found SAF for: $GROUP_NAME"
            fi
        done
    done
    echo "Found ${#ANGSD_JOBS[@]} groups with existing SAF files"
fi

# -------------------------------
# Step 2c: PCAngsd
# -------------------------------
if [[ $RUN_PCA -eq 1 ]]; then
    echo ""
    echo "==================================================================="
    echo "STEP 2c: Running PCAngsd"
    echo "==================================================================="

    if [[ -n "${LD_PRUNE_JOB:-}" ]]; then
        PCA_DEP="--dependency=afterok:$LD_PRUNE_JOB"
    else
        PCA_DEP=""
    fi

    PCANGSD_JOB=$(sbatch --parsable \
        $PCA_DEP \
        --job-name="pcangsd_ALL" \
        --cpus-per-task=40 \
        --mem=160G \
        --time=96:00:00 \
        --export=ALL,BAM_LIST="$ALL_LIST",REFERENCE="$REFERENCE",OUT_DIR="$OUT_DIR",THREADS=40 \
        scripts/p06c3_pcangsd.sh)

    echo "  Submitted: PCAngsd Job $PCANGSD_JOB"
fi


# Step 2d: STRUCTURE preparation
if [[ $RUN_PCA -eq 1 ]]; then
    echo ""
    echo "============================"
    echo "STEP 2d: STRUCTURE prep"
    echo "============================"
    
    STRUCT_JOB=$(sbatch --parsable \
        --dependency=afterok:$PCA_JOB \
        --job-name="structure_prep" \
        --export=ALL,OUT_DIR="$OUT_DIR",\
BEAGLE_FILE="$OUT_DIR/pop_map_ALL.beagle.gz",\
METADATA="$OUT_DIR/metadata.tsv" \
        scripts/p06g_structure.sh)
    
    echo "  Submitted: STRUCTURE prep $STRUCT_JOB"
fi


# -------------------------------
# Step 3: Theta calculation
# -------------------------------
if [[ $RESUME_FROM -le 3 ]]; then
    echo ""
    echo "==================================================================="
    echo "STEP 3: Calculating theta per group"
    echo "==================================================================="

    declare -A THETA_JOBS

    for GROUP_NAME in "${!ANGSD_JOBS[@]}"; do
        DEP_JOB=${ANGSD_JOBS[$GROUP_NAME]}
        
        LIST_FILE=""
        for LIST_DIR in "pop_map_lists" "site_year_map_lists" "year_map_lists"; do
            CANDIDATE="$OUT_DIR/$LIST_DIR/${GROUP_NAME}.list"
            [[ -f "$CANDIDATE" ]] && LIST_FILE="$CANDIDATE" && break
        done
        
        [[ -z "$LIST_FILE" ]] && echo "WARNING: List file missing for $GROUP_NAME, skipping theta" && continue
        
        # Skip dependency if resuming (DEP_JOB="completed")
        if [[ "$DEP_JOB" == "completed" ]]; then
            JOB_ID=$(sbatch --parsable \
                --job-name="theta_${GROUP_NAME}" \
                --cpus-per-task=4 \
                --mem=16G \
                --time=06:00:00 \
                --export=ALL,BAM_LIST="$LIST_FILE",REFERENCE="$REFERENCE",OUT_DIR="$OUT_DIR",ANCESTRAL="$ANCESTRAL",MIN_MAPQ="$MIN_MAPQ",MIN_Q="$MIN_Q",FOLD_SFS="$FOLD_SFS",WINDOW_SIZE="$WINDOW_SIZE",WINDOW_STEP="$WINDOW_STEP" \
                scripts/p06d_merge_theta.sh)
        else
            JOB_ID=$(sbatch --parsable \
                --dependency=afterok:$DEP_JOB \
                --job-name="theta_${GROUP_NAME}" \
                --cpus-per-task=4 \
                --mem=16G \
                --time=06:00:00 \
                --export=ALL,BAM_LIST="$LIST_FILE",REFERENCE="$REFERENCE",OUT_DIR="$OUT_DIR",ANCESTRAL="$ANCESTRAL",MIN_MAPQ="$MIN_MAPQ",MIN_Q="$MIN_Q",FOLD_SFS="$FOLD_SFS",WINDOW_SIZE="$WINDOW_SIZE",WINDOW_STEP="$WINDOW_STEP" \
                scripts/p06d_merge_theta.sh)
        fi
        
        THETA_JOBS["$GROUP_NAME"]=$JOB_ID
    done

    echo "Total theta jobs: ${#THETA_JOBS[@]}"
else
    echo "Skipping Step 3 (theta calculation)"
    declare -A THETA_JOBS
fi

# -------------------------------
# Step 4: Pairwise FST
# -------------------------------
if [[ $RESUME_FROM -le 4 ]]; then
    echo ""
    echo "==================================================================="
    echo "STEP 4: Calculating pairwise FST"
    echo "==================================================================="

    FST_JOBS=()

    submit_fst() {
        local L1=$1; local L2=$2
        local P1=$(basename "$L1" .list)
        local P2=$(basename "$L2" .list)
        
        THETA1="$OUT_DIR/${P1}.thetas.idx"
        THETA2="$OUT_DIR/${P2}.thetas.idx"

        if [[ ! -f "$THETA1" || ! -f "$THETA2" ]]; then
            echo "Skipping $P1 vs $P2 (theta files missing)"
            return
        fi
        
        N1=$(wc -l < "$L1"); N2=$(wc -l < "$L2")
        ((N1<2||N2<2)) && return
        
        if [[ -n "${THETA_JOBS[$P1]:-}" && -n "${THETA_JOBS[$P2]:-}" ]]; then
            DEP_JOBS="${THETA_JOBS[$P1]}:${THETA_JOBS[$P2]}"
            DEP="--dependency=afterok:$DEP_JOBS"
        else
            DEP=""
        fi
        
        JOB_ID=$(sbatch --parsable \
            $DEP \
            --job-name="fst_${P1}_${P2}" \
            --cpus-per-task=8 \
            --mem=228G \
            --time=12:00:00 \
            --export=ALL,LIST1="$L1",LIST2="$L2",REFERENCE="$REFERENCE",OUT_DIR="$OUT_DIR",ANCESTRAL="$ANCESTRAL",THREADS=8,MIN_MAPQ="$MIN_MAPQ",MIN_Q="$MIN_Q",FOLD_SFS="$FOLD_SFS",MIN_IND_RATIO="$MIN_IND_RATIO",WINDOW_SIZE="$WINDOW_SIZE",WINDOW_STEP="$WINDOW_STEP" \
            scripts/p06e_fst.sh)
        
        echo "    $P1 vs $P2: Job $JOB_ID"
        FST_JOBS+=("$JOB_ID")
    }

    for DIR in "pop_map_lists" "site_year_map_lists" "year_map_lists"; do
        [[ ! -d "$OUT_DIR/$DIR" ]] && continue
        LISTS=("$OUT_DIR/$DIR"/*.list)
        for ((i=0;i<${#LISTS[@]}-1;i++)); do
            for ((j=i+1;j<${#LISTS[@]};j++)); do
                submit_fst "${LISTS[i]}" "${LISTS[j]}"
            done
        done
    done

    echo ""
    echo "Total FST jobs: ${#FST_JOBS[@]}"
else
    echo "Skipping Step 4 (FST calculation)"
    FST_JOBS=()
fi

# -------------------------------
# Step 5: Summary
# -------------------------------
if [[ $RESUME_FROM -le 5 && ${#FST_JOBS[@]} -gt 0 ]]; then
    echo ""
    echo "==================================================================="
    echo "STEP 5: Generating summary statistics"
    echo "==================================================================="

    DEP_JOBS=$(IFS=:; echo "${FST_JOBS[*]}")
    
    SUMMARY_JOB=$(sbatch --parsable \
        --dependency=afterok:$DEP_JOBS \
        --job-name="summary" \
        --cpus-per-task=1 \
        --mem=4G \
        --time=01:00:00 \
        --export=ALL,OUT_DIR="$OUT_DIR",WINDOW_SIZE="$WINDOW_SIZE",WINDOW_STEP="$WINDOW_STEP" \
        scripts/p06f_summary.sh)
    
    echo "  Submitted: Job $SUMMARY_JOB"
fi

echo ""
echo "==================================================================="
echo "Pipeline Submission Complete"
echo "==================================================================="
echo ""
echo "Job Summary:"
echo "  Step 0 (Subsampled ANGSD for LD):  ${ALL_ANGSD_JOB:-skipped}"
echo "  Step 1 (LD pruning):                ${LD_PRUNE_JOB:-skipped}"
echo "  Step 2 (SAF per group):             ${NUM_GROUPS:-0} jobs"
echo "  Step 2b (PCA for ALL):              ${PCA_JOB:-skipped (use -P 1)}"
echo "  Step 3 (Theta):                     ${#THETA_JOBS[@]} jobs"
echo "  Step 4 (FST):                       ${#FST_JOBS[@]} jobs"
echo "  Step 5 (Summary):                   ${SUMMARY_JOB:-pending}"
echo ""
echo "Output structure:"
echo "  $OUT_DIR/metadata.tsv           - Sample metadata"
echo "  $OUT_DIR/pop_map_lists/         - By population + ALL"
echo "  $OUT_DIR/site_year_map_lists/   - By site-year"
echo "  $OUT_DIR/year_map_lists/        - By year"
echo "  $OUT_DIR/LDpruned_snps.list     - LD-pruned SNP list"
echo ""
echo "Monitor:"
echo "  squeue -u \$USER"
echo "  tail -f $OUT_DIR/logs/*.log"
echo "==================================================================="