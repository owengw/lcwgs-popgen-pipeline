#!/bin/bash
#SBATCH --job-name=p01_trim_array
#SBATCH --output=p01_trim_array_%A_%a.log
#SBATCH --mem=40G
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --time=48:00:00
#SBATCH --array=1-333%20   # Adjust to number of rows in your prefix list

source ~/.bash_profile
conda activate /mnt/parscratch/users/bi4og/conda_envs/owenspopgen

helpFunction() {
    echo ""
    echo "Usage: $0 -p <raw_data_path> -l <prefix_list> -f <forward_ext> -r <reverse_ext> -k <ILLUMINACLIP> -s <SLIDINGWINDOW> -L <LEADING> -t <TRAILING> -m <MINLEN>"
    echo ""
    echo "  -p Path to raw FASTQ files (default: raw_data)"
    echo "  -l File containing prefixes, one per line (default: physalia/prefix_list.txt)"
    echo "  -f Forward reads extension (default: _R1.fastq.gz)"
    echo "  -r Reverse reads extension (default: _R2.fastq.gz)"
    echo "  -k ILLUMINACLIP parameters (default: ILLUMINACLIP:scripts/NEBNext.fa:2:30:10)"
    echo "  -s SLIDINGWINDOW parameters (default: SLIDINGWINDOW:4:20)"
    echo "  -L LEADING parameter (default: LEADING:3)"
    echo "  -t TRAILING parameter (default: TRAILING:3)"
    echo "  -m MINLEN parameter (default: MINLEN:36)"
    exit 1
}

# Parse options
while getopts "p:l:f:r:k:s:L:t:m:h" opt; do
    case "$opt" in
        p) RAWPATH="$OPTARG" ;;
        l) PREFIXLIST="$OPTARG" ;;
        f) parameterF="$OPTARG" ;;
        r) parameterR="$OPTARG" ;;
        k) parameterK="$OPTARG" ;;
        s) parameterS="$OPTARG" ;;
        L) parameterL="$OPTARG" ;;
        t) parameterT="$OPTARG" ;;
        m) parameterM="$OPTARG" ;;
        h) helpFunction ;;
        *) helpFunction ;;
    esac
done

# Defaults
RAWPATH=${RAWPATH:-raw_data}
PREFIXLIST=${PREFIXLIST:-physalia/prefix_list.txt}
parameterF=${parameterF:-_R1.fastq.gz}
parameterR=${parameterR:-_R2.fastq.gz}
parameterK=${parameterK:-ILLUMINACLIP:scripts/NEBNext.fa:2:30:10}
parameterS=${parameterS:-SLIDINGWINDOW:4:20}
parameterL=${parameterL:-LEADING:3}
parameterT=${parameterT:-TRAILING:3}
parameterM=${parameterM:-MINLEN:36}

src=$PWD
mkdir -p "$src/physalia/adapter_clipped"

# Get the prefix for this array task
FILE_PREFIX=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$PREFIXLIST")

if [ -z "$FILE_PREFIX" ]; then
    echo "ERROR: No prefix found for array task $SLURM_ARRAY_TASK_ID"
    exit 1
fi

echo "=========================================="
echo "Processing prefix: $FILE_PREFIX"
echo "Array ID:        $SLURM_ARRAY_TASK_ID"
echo "Time:            $(date)"
echo "=========================================="
echo ""

# Build input paths from the prefix
R1_FILE="$src/$RAWPATH/${FILE_PREFIX}${parameterF}"
R2_FILE="$src/$RAWPATH/${FILE_PREFIX}${parameterR}"

# Determine if paired-end or single-end from existence of R2
if [[ -f "$R2_FILE" ]]; then
    DATATYPE="pe"
else
    DATATYPE="se"
fi

# Output base
SAMPLE_UNIQ_ID=$(basename "$FILE_PREFIX")
SAMPLEADAPT="$src/physalia/adapter_clipped/${SAMPLE_UNIQ_ID}"
mkdir -p "$(dirname "$SAMPLEADAPT")"

# Trimmomatic
TRIMMOMATIC_JAR="$CONDA_PREFIX/share/trimmomatic-0.40-0/trimmomatic.jar"
JAVA_OPTS="-Xmx10G"

if [[ "$DATATYPE" == "pe" ]]; then
    java $JAVA_OPTS -jar "$TRIMMOMATIC_JAR" PE \
        -threads ${SLURM_CPUS_PER_TASK} -phred33 \
        "$R1_FILE" "$R2_FILE" \
        "${SAMPLEADAPT}_f_paired.fastq.gz" \
        "${SAMPLEADAPT}_f_unpaired.fastq.gz" \
        "${SAMPLEADAPT}_r_paired.fastq.gz" \
        "${SAMPLEADAPT}_r_unpaired.fastq.gz" \
        "$parameterK" "$parameterS" "$parameterL" "$parameterT" "$parameterM"
else
    java $JAVA_OPTS -jar "$TRIMMOMATIC_JAR" SE \
        -threads ${SLURM_CPUS_PER_TASK} -phred33 \
        "$R1_FILE" \
        "${SAMPLEADAPT}_se.fastq.gz" \
        "$parameterK" "$parameterS" "$parameterL" "$parameterT" "$parameterM"
fi

echo ""
echo "=========================================="
echo "COMPLETE: $SAMPLE_UNIQ_ID"
echo "Time: $(date)"
echo "=========================================="
