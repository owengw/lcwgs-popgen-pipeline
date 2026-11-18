#!/bin/bash
#SBATCH --job-name=p04_merge_bam_array
#SBATCH --output=p04_merge_bam_array_%A_%a.log
#SBATCH --mem=32G
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --time=24:00:00
#SBATCH --array=1-198%20  # number of unique biological samples

source ~/.bash_profile
conda activate /mnt/parscratch/users/bi4og/conda_envs/owenspopgen

helpFunction() {
    echo ""
    echo "Usage: $0 -q <prefix_list> -i <input_dir> -o <output_dir> -r <reference_name>"
    echo ""
    echo "  -q  File with full prefixes including lanes (default: physalia/prefix_list.txt)"
    echo "  -i  Input BAM directory (default: physalia/aligned)"
    echo "  -o  Output BAM directory (default: physalia/merged)"
    echo "  -r  Reference name used in BAM filenames (default: Plong_genome_flye)"
    echo ""
    exit 1
}

# Parse options
while getopts "q:i:o:r:h" opt; do
    case "$opt" in
        q) parameterQ="$OPTARG" ;;
        i) parameterI="$OPTARG" ;;
        o) parameterO="$OPTARG" ;;
        r) parameterR="$OPTARG" ;;
        h) helpFunction ;;
        *) helpFunction ;;
    esac
done

# Defaults
parameterQ=${parameterQ:-physalia/prefix_list.txt}
parameterI=${parameterI:-physalia/aligned}
parameterO=${parameterO:-physalia/merged}
parameterR=${parameterR:-Plong_genome_flye}

src=$PWD
INPUT_DIR="$src/$parameterI"
OUTPUT_DIR="$src/$parameterO"
REFERENCE="$parameterR"

BAMLIST_DIR="$src/physalia/sample_lists"
BAMLIST_FILE="$BAMLIST_DIR/bam_list_dedup_overlapclipped.list"
LOCKFILE="$BAMLIST_FILE.lock"

# Calculate Picard memory (leave 4G headroom)
TOTAL_MEM_MB=$((${SLURM_MEM_PER_NODE:-32768}))
PICARD_MEM_MB=$((TOTAL_MEM_MB - 4096))
JAVA_OPTS="-Xmx${PICARD_MEM_MB}m"

mkdir -p "$OUTPUT_DIR" "$BAMLIST_DIR" "$OUTPUT_DIR/tmp"

# Extract unique biological sample IDs
# Everything before first underscore
# Example:
#   282-NS7_CTGAACGTAT-GCCAATACAT_L007 -> 282-NS7
UNIQUE_SAMPLES=($(awk -F_ '{print $1}' "$src/$parameterQ" | sort -u))
SAMPLE_TAG="${UNIQUE_SAMPLES[$((SLURM_ARRAY_TASK_ID-1))]}"

if [ -z "$SAMPLE_TAG" ]; then
    echo "ERROR: No sample found for array task ${SLURM_ARRAY_TASK_ID}"
    exit 1
fi

echo "=========================================="
echo "Processing biological sample: ${SAMPLE_TAG}"
echo "Array task ID:               ${SLURM_ARRAY_TASK_ID}"
echo "Time:                        $(date)"
echo "=========================================="
echo ""

# Find BAMs across all lanes for this sample
PE_BAMS=($(ls "$INPUT_DIR/${SAMPLE_TAG}"*_pe_bt2_${REFERENCE}_minq20_sorted.bam 2>/dev/null))
SE_BAMS=($(ls "$INPUT_DIR/${SAMPLE_TAG}"*_se_bt2_${REFERENCE}_minq20_sorted.bam 2>/dev/null))

echo "Found ${#PE_BAMS[@]} PE BAM(s) and ${#SE_BAMS[@]} SE BAM(s)"


# Function to merge, dedup, clip, index, and append to BAM list
process_bams() {
    local TYPE="$1"
    shift
    local BAM_ARRAY=("$@")

    local NUM_FILES=${#BAM_ARRAY[@]}
    if [ "$NUM_FILES" -eq 0 ]; then
        echo "No $TYPE BAMs found for sample $SAMPLE_TAG"
        return
    fi

    echo ""
    echo "Processing $NUM_FILES $TYPE BAM(s) for $SAMPLE_TAG"
    for f in "${BAM_ARRAY[@]}"; do
        echo "  [FOUND] $(basename "$f")"
    done
    echo ""

    local MERGED_BAM="${OUTPUT_DIR}/${SAMPLE_TAG}_merged_${TYPE}_bt2_${REFERENCE}_minq20_sorted.bam"

    if [ "$NUM_FILES" -eq 1 ]; then
        cp "${BAM_ARRAY[0]}" "$MERGED_BAM"
    else
        samtools merge -@ ${SLURM_CPUS_PER_TASK} "$MERGED_BAM" "${BAM_ARRAY[@]}"
    fi
    [ -f "$MERGED_BAM" ] || { echo "ERROR: merge failed"; exit 1; }

    # Deduplicate
    local DEDUP_BAM="${OUTPUT_DIR}/${SAMPLE_TAG}_merged_${TYPE}_bt2_${REFERENCE}_minq20_sorted_dedup.bam"
    local DUPSTAT="${OUTPUT_DIR}/${SAMPLE_TAG}_merged_${TYPE}_bt2_${REFERENCE}_minq20_sorted_dupstat.txt"

    picard $JAVA_OPTS MarkDuplicates \
        I="$MERGED_BAM" \
        O="$DEDUP_BAM" \
        M="$DUPSTAT" \
        REMOVE_DUPLICATES=true \
        MAX_RECORDS_IN_RAM=1000000 \
        CREATE_INDEX=false \
        CREATE_MD5_FILE=false \
        VALIDATION_STRINGENCY=LENIENT \
        TMP_DIR="$OUTPUT_DIR/tmp"

    [ -f "$DEDUP_BAM" ] || { echo "ERROR: MarkDuplicates failed"; exit 1; }

    # Clip overlaps
    local FINAL_BAM="${OUTPUT_DIR}/${SAMPLE_TAG}_merged_${TYPE}_bt2_${REFERENCE}_minq20_sorted_dedup_overlapclipped.bam"

    bam clipOverlap \
        --in "$DEDUP_BAM" \
        --out "$FINAL_BAM" \
        --stats

    [ -f "$FINAL_BAM" ] || { echo "ERROR: clipOverlap failed"; exit 1; }

    # Index + sanity check
    samtools index "$FINAL_BAM"
    samtools quickcheck "$FINAL_BAM" || { echo "ERROR: samtools quickcheck failed"; exit 1; }
    [ -s "$FINAL_BAM" ] || { echo "ERROR: BAM is empty"; exit 1; }

    # Append to BAM list (array-safe)
    (
        flock -x 200
        echo "$FINAL_BAM" >> "$BAMLIST_FILE"
    ) 200>"$LOCKFILE"

    # Cleanup intermediates
    rm -f "$MERGED_BAM" "$DEDUP_BAM"

    echo ""
    echo "=========================================="
    echo "COMPLETE: ${SAMPLE_TAG} ($TYPE)"
    echo "Final BAM: $(basename "$FINAL_BAM")"
    du -h "$FINAL_BAM" | awk '{print "Size: " $1}'
    echo "Time: $(date)"
    echo "=========================================="
}


# Run PE and SE separately
process_bams "pe" "${PE_BAMS[@]}"
process_bams "se" "${SE_BAMS[@]}"
