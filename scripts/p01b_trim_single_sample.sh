#!/bin/bash
set -euo pipefail

TRIMMOMATIC_JAR="/mnt/parscratch/users/bi4og/conda_envs/owenspopgen/share/trimmomatic-0.40-0/trimmomatic.jar"
JAVA_OPTS="-Xmx10G" #sets the java memory to 8G so it doesn't run out

# Positional arguments from parallel
FILE_PREFIX="$1"
PARAMETERF="$2"
PARAMETERR="$3"
PARAMETERK="$4"
PARAMETERS="$5"
PARAMETERL="$6"
PARAMETERT="$7"
PARAMETERM="${8}"
PARAMETERQ="${9}"
PARAMETERP="${10}"

src=$PWD

TABLE_PREFIX="${FILE_PREFIX%%-*}"

# Debug
echo "DEBUG: FILE_PREFIX=$FILE_PREFIX" >&2
echo "DEBUG: Full path Q=$src/$PARAMETERQ" >&2
echo "DEBUG: Full path P=$src/$PARAMETERP" >&2

# Extract sample info
read SAMPLE_ID POP_ID SEQ_ID LANE_ID DATATYPE < <(
    awk -F'\t' -v p="$TABLE_PREFIX" '$1==p {print $4,$5,$3,$2,$6; exit}' "$src/$PARAMETERQ"
)

SAMPLE_UNIQ_ID="${SAMPLE_ID}_${POP_ID}_${SEQ_ID}_${LANE_ID}"
SAMPLEADAPT="$src/physalia/adapter_clipped/$SAMPLE_UNIQ_ID"
mkdir -p "$(dirname "$SAMPLEADAPT")"

# Run Trimmomatic safely with quotes
if [[ "$DATATYPE" == "pe" ]]; then
    java $JAVA_OPTS -jar "$TRIMMOMATIC_JAR" PE \
        -threads 4 -phred33 \
        "$src/$PARAMETERP/${FILE_PREFIX}${PARAMETERF}" \
        "$src/$PARAMETERP/${FILE_PREFIX}${PARAMETERR}" \
        "${SAMPLEADAPT}_adapter_clipped_f_paired.fastq.gz" \
        "${SAMPLEADAPT}_adapter_clipped_f_unpaired.fastq.gz" \
        "${SAMPLEADAPT}_adapter_clipped_r_paired.fastq.gz" \
        "${SAMPLEADAPT}_adapter_clipped_r_unpaired.fastq.gz" \
        "$PARAMETERK" "$PARAMETERS" "$PARAMETERL" "$PARAMETERT" "$PARAMETERM"
        #-trimlog "${SAMPLEADAPT}_trimlog.txt"

elif [[ "$DATATYPE" == "se" ]]; then
    java $JAVA_OPTS -jar "$TRIMMOMATIC_JAR" SE \
        -threads 4 -phred33 \
        "$src/$PARAMETERP/${FILE_PREFIX}${PARAMETERF}" \
        "${SAMPLEADAPT}_se.fastq.gz" \
        "$PARAMETERK" "$PARAMETERS" "$PARAMETERL" "$PARAMETERT" "$PARAMETERM"
        #-trimlog "${SAMPLEADAPT}_trimlog.txt"
else
    echo "ERROR: Invalid DATATYPE '$DATATYPE'" >&2
    exit 1
fi
