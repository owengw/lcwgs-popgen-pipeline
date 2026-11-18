#!/bin/bash
#SBATCH --job-name=p05_depth
#SBATCH --output=p05_depth_%A_%a.log
#SBATCH --mem=16G
#SBATCH --cpus-per-task=4
#SBATCH --time=12:00:00
#SBATCH --array=1-198%20

source ~/.bash_profile
conda activate /mnt/parscratch/users/bi4og/conda_envs/owenspopgen

helpFunction() {
    echo ""
    echo "Usage: $0 -l <bam_list> -o <output_dir> -q <min_mapq> -Q <min_baseq>"
    echo ""
    echo "  -l  BAM list file (deduplicated, overlap-clipped BAMs)"
    echo "  -o  Output directory for depth files"
    echo "  -q  Minimum mapping quality"
    echo "  -Q  Minimum base quality"
    echo ""
    exit 1
}

# Parse options
while getopts "l:o:q:Q:h" opt; do
    case "$opt" in
        l) BAM_LIST="$OPTARG" ;;
        o) DEPTH_DIR="$OPTARG" ;;
        q) MIN_MAPQ="$OPTARG" ;;
        Q) MIN_BASEQ="$OPTARG" ;;
        h) helpFunction ;;
        *) helpFunction ;;
    esac
done

# Defaults - updated to point to deduplicated BAMs (no indel realignment)
src=$PWD
BAM_LIST=${BAM_LIST:-$src/physalia/sample_lists/bam_list_dedup_overlapclipped.list}
DEPTH_DIR=${DEPTH_DIR:-$src/physalia/depths}
MIN_MAPQ=${MIN_MAPQ:-20}
MIN_BASEQ=${MIN_BASEQ:-20}

mkdir -p "$DEPTH_DIR"/{depth_files,stats}

# Fetch BAM for this array task
BAM_FILE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$BAM_LIST")
[ -n "$BAM_FILE" ] || { echo "ERROR: No BAM for array task $SLURM_ARRAY_TASK_ID"; exit 1; }

# Handle both absolute and relative paths
[[ "$BAM_FILE" = /* ]] || BAM_FILE="$src/$BAM_FILE"
[ -f "$BAM_FILE" ] || { echo "ERROR: BAM not found: $BAM_FILE"; exit 1; }

SAMPLE=$(basename "$BAM_FILE" _dedup_overlapclipped.bam)
DEPTH_FILE="$DEPTH_DIR/depth_files/${SAMPLE}.depth.gz"
STATS_FILE="$DEPTH_DIR/stats/${SAMPLE}.depth_stats.tsv"

echo "=========================================="
echo "Sample:        $SAMPLE"
echo "BAM:           $BAM_FILE"
echo "Depth file:    $DEPTH_FILE"
echo "Stats file:    $STATS_FILE"
echo "MAPQ >=        $MIN_MAPQ"
echo "BaseQ >=       $MIN_BASEQ"
echo "Array task:    $SLURM_ARRAY_TASK_ID"
echo "Time:          $(date)"
echo "=========================================="
echo ""

# Depth calculation
echo "Calculating depth..."
samtools depth -aa -q "$MIN_MAPQ" -Q "$MIN_BASEQ" -@ "$SLURM_CPUS_PER_TASK" "$BAM_FILE" | \
gzip > "$DEPTH_FILE"

if [ ! -f "$DEPTH_FILE" ]; then
    echo "ERROR: Depth file not created"
    exit 1
fi

# Check if depth file is not empty
if [ ! -s "$DEPTH_FILE" ]; then
    echo "ERROR: Depth file is empty"
    exit 1
fi

echo "Depth file created successfully"

# Per-sample stats
echo "Calculating statistics..."
zcat "$DEPTH_FILE" | \
awk -v S="$SAMPLE" '
{
    d=$3; n++; sum+=d; sumsq+=d*d
    if(d>0) cov++
    if(d>=1) x1++; if(d>=5) x5++; if(d>=10) x10++; if(d>=20) x20++; if(d>=30) x30++
    if(d>max) max=d
}
END {
    if(n==0) {
        print "ERROR: No data processed" > "/dev/stderr"
        exit 1
    }
    mean=sum/n
    sd=sqrt(sumsq/n - mean^2)
    prop_cov = (n>0) ? cov/n : 0
    printf "sample\ttotal_sites\tmean_depth\tsd_depth\tmax_depth\tsites_covered\tprop_covered\tsites_1x\tsites_5x\tsites_10x\tsites_20x\tsites_30x\n"
    printf "%s\t%d\t%.3f\t%.3f\t%d\t%d\t%.4f\t%d\t%d\t%d\t%d\t%d\n", S,n,mean,sd,max,cov,prop_cov,x1,x5,x10,x20,x30
}' > "$STATS_FILE"

if [ ! -s "$STATS_FILE" ]; then
    echo "ERROR: Stats file is empty or failed"
    exit 1
fi

echo "Sample complete: $SAMPLE"
echo "=========================================="