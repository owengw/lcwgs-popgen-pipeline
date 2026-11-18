#!/bin/bash
#SBATCH --job-name=p_sex_identification
#SBATCH --output=p_sex_identification_%a_%j.log
#SBATCH --mem=16G
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=4
#SBATCH --time=04:00:00
#SBATCH --array=1-20

source ~/.bash_profile
conda activate /mnt/parscratch/users/bi4og/conda_envs/owenspopgen

src=/mnt/parscratch/users/bi4og
BAMDIR=$src/physalia/merged
OUTDIR=$src/sex_identification
AGP=$src/ragtag_output/ragtag.scaffold.agp

mkdir -p $OUTDIR

# ============================================================
# STEP 1: Generate contig lists and BED files (task 1 only)
# ============================================================
if [ "$SLURM_ARRAY_TASK_ID" -eq 1 ]; then
    echo "[$(date)] Generating contig lists from AGP..."

    # Chr5 contigs (putative X)
    awk '$1=="chr5_RagTag" && $5=="W" {print $6}' $AGP \
        > $OUTDIR/chr5_contigs.txt

    # Autosomal contigs (chr1-4, chr6-13)
    awk '$5=="W" && $1!~/chr5_RagTag/ && $1~/^chr[0-9]*_RagTag$/' $AGP | \
        awk '{print $6}' \
        > $OUTDIR/autosomal_contigs.txt

    echo "Chr5 contigs:      $(wc -l < $OUTDIR/chr5_contigs.txt)"
    echo "Autosomal contigs: $(wc -l < $OUTDIR/autosomal_contigs.txt)"

    # Build BED files from BAM header using contig lists
    EXAMPLE_BAM=$(ls $BAMDIR/*_overlapclipped.bam | head -1)

    # Replace the BED generation lines with:
    samtools view -H $EXAMPLE_BAM | grep "^@SQ" | \
        sed 's/@SQ\tSN:\([^\t]*\)\tLN:\([0-9]*\)/\1\t0\t\2/' | \
        awk 'NR==FNR{contigs[$1]=1; next} $1 in contigs' \
        $OUTDIR/chr5_contigs.txt - \
        > $OUTDIR/chr5_regions.bed

    samtools view -H $EXAMPLE_BAM | grep "^@SQ" | \
        sed 's/@SQ\tSN:\([^\t]*\)\tLN:\([0-9]*\)/\1\t0\t\2/' | \
        awk 'NR==FNR{contigs[$1]=1; next} $1 in contigs' \
        $OUTDIR/autosomal_contigs.txt - \
        > $OUTDIR/autosomal_regions.bed

    echo "Chr5 BED regions:      $(wc -l < $OUTDIR/chr5_regions.bed)"
    echo "Autosomal BED regions: $(wc -l < $OUTDIR/autosomal_regions.bed)"

    # Verify BED files have content
    if [ ! -s $OUTDIR/chr5_regions.bed ]; then
        echo "ERROR: chr5_regions.bed is empty - contig name mismatch?"
        head -3 $OUTDIR/chr5_contigs.txt
        samtools view -H $EXAMPLE_BAM | grep "^@SQ" | head -3
        exit 1
    fi
else
    # Wait for task 1 to generate BED files
    echo "[$(date)] Task $SLURM_ARRAY_TASK_ID waiting for BED files..."
    for i in $(seq 1 12); do
        [ -f $OUTDIR/chr5_regions.bed ] && break
        sleep 10
    done
    if [ ! -f $OUTDIR/chr5_regions.bed ]; then
        echo "ERROR: BED files not generated after 120s"
        exit 1
    fi
fi

# ============================================================
# STEP 2: Calculate coverage per individual
# ============================================================
mapfile -t BAMS < <(ls $BAMDIR/*_overlapclipped.bam | sort)
TOTAL=${#BAMS[@]}
BATCH_SIZE=$(( (TOTAL + 19) / 20 ))
START=$(( (SLURM_ARRAY_TASK_ID - 1) * BATCH_SIZE ))
END=$(( START + BATCH_SIZE - 1 ))

echo "[$(date)] Task $SLURM_ARRAY_TASK_ID: indices $START-$END of $TOTAL"

TASK_OUT=$OUTDIR/sex_coverage_task${SLURM_ARRAY_TASK_ID}.csv

for idx in $(seq $START $END); do
    [ $idx -ge $TOTAL ] && break
    BAM=${BAMS[$idx]}
    SAMPLE=$(basename $BAM \
        _merged_pe_bt2_Plong_genome_flye_minq20_sorted_dedup_overlapclipped.bam)

    echo "[$(date)] Processing: $SAMPLE"

    # Mean depth on chr5 contigs
    CHR5_COV=$(samtools depth -a \
        -b $OUTDIR/chr5_regions.bed \
        $BAM 2>/dev/null | \
        awk '{sum+=$3; n++} END{if(n>0) printf "%.4f", sum/n; else print 0}')

    # Mean depth on autosomal contigs
    AUTO_COV=$(samtools depth -a \
        -b $OUTDIR/autosomal_regions.bed \
        $BAM 2>/dev/null | \
        awk '{sum+=$3; n++} END{if(n>0) printf "%.4f", sum/n; else print 0}')

    # X/autosome ratio
    RATIO=$(awk "BEGIN{
        if ($AUTO_COV > 0)
            printf \"%.4f\", $CHR5_COV / $AUTO_COV
        else
            print \"NA\"
    }")

    echo "$SAMPLE,$CHR5_COV,$AUTO_COV,$RATIO" >> $TASK_OUT
    echo "[$(date)] $SAMPLE: chrX=$CHR5_COV auto=$AUTO_COV ratio=$RATIO"
done

echo "[$(date)] Task $SLURM_ARRAY_TASK_ID complete"

# ============================================================
# STEP 3: Merge all task outputs (task 20 only)
# ============================================================
if [ "$SLURM_ARRAY_TASK_ID" -eq 20 ]; then
    echo "[$(date)] Waiting 120s for other tasks..."
    sleep 120

    echo "sample,chrX_cov,auto_cov,ratio" > $OUTDIR/sex_coverage.csv
    cat $OUTDIR/sex_coverage_task*.csv | sort >> $OUTDIR/sex_coverage.csv

    echo "[$(date)] Merged: $OUTDIR/sex_coverage.csv"
    echo "Total individuals: $(tail -n +2 $OUTDIR/sex_coverage.csv | wc -l)"

    # Quick sanity check
    echo "Coverage summary:"
    awk -F',' 'NR>1 && $4!="NA" {sum+=$4; n++}
               END{printf "Mean ratio: %.3f (n=%d)\n", sum/n, n}' \
        $OUTDIR/sex_coverage.csv
fi