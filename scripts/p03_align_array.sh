#!/bin/bash
#SBATCH --job-name=p03_align_array
#SBATCH --output=logs/p03_align_array_%A_%a.log
#SBATCH --mem=16G
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --time=24:00:00
#SBATCH --array=1-333%20     # set to number of lines in your prefix list

# FIXED: removed broken line-continuations (trailing spaces after '\') and the
# split samtools view|sort pipe that made the original fail. Default genome now
# the annotated RagTag/pilon assembly. set -eu (no pipefail).
# Usage: sbatch scripts/p03_align_array.sh -g <genome> -l <prefix_list> [-p <adapter_clipped>]

source ~/.bash_profile
conda activate /mnt/parscratch/users/bi4og/conda_envs/owenspopgen
set -eu

helpFunction() {
  echo "Usage: $0 -g <genome_fasta> -l <prefix_list> -p <adapter_clipped_path>"
  exit 1
}
while getopts "g:l:p:h" opt; do
  case "$opt" in
    g) GENOME="$OPTARG" ;;
    l) PREFIXLIST="$OPTARG" ;;
    p) ADAPT_CLIP_PATH="$OPTARG" ;;
    h) helpFunction ;;
    *) helpFunction ;;
  esac
done

GENOME=${GENOME:-ncbi_submission/plongirostris.fsa}
PREFIXLIST=${PREFIXLIST:-physalia/prefix_list.txt}
ADAPT_CLIP_PATH=${ADAPT_CLIP_PATH:-physalia/adapter_clipped}
src=$PWD
mkdir -p "$src/physalia/aligned" logs

REFERENCE="$src/$GENOME"
REFBASENAME="${REFERENCE%.*}"
REFNAME=$(basename "$REFBASENAME")

PREFIX=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$PREFIXLIST")
[ -n "$PREFIX" ] || { echo "ERROR: no prefix for array task $SLURM_ARRAY_TASK_ID"; exit 1; }

SAMPLE_UNIQ_ID=$(basename "$PREFIX")
SAMPLETOMAP="$src/$ADAPT_CLIP_PATH/$SAMPLE_UNIQ_ID"
SAMPLEBAM="$src/physalia/aligned/$SAMPLE_UNIQ_ID"
FASTQSUFFIX1=_f_paired.fastq.gz
FASTQSUFFIX2=_r_paired.fastq.gz

if [[ -f "${SAMPLETOMAP}${FASTQSUFFIX2}" ]]; then DATATYPE="pe"; else DATATYPE="se"; fi
echo "Processing: $SAMPLE_UNIQ_ID | datatype: $DATATYPE | ref: $REFNAME | $(date)"
PU=$SAMPLE_UNIQ_ID

if [[ "$DATATYPE" == "pe" ]]; then
  for f in "${SAMPLETOMAP}${FASTQSUFFIX1}" "${SAMPLETOMAP}${FASTQSUFFIX2}"; do
    [ -f "$f" ] || { echo "ERROR: file not found: $f"; exit 1; }
  done
else
  [ -f "${SAMPLETOMAP}${FASTQSUFFIX1}" ] || { echo "ERROR: file not found: ${SAMPLETOMAP}${FASTQSUFFIX1}"; exit 1; }
fi

MAPPINGPRESET=very-sensitive
OUTSAM="${SAMPLEBAM}_${DATATYPE}_bt2_${REFNAME}.sam"
OUTBAM="${SAMPLEBAM}_${DATATYPE}_bt2_${REFNAME}.bam"
SORTED="${SAMPLEBAM}_${DATATYPE}_bt2_${REFNAME}_minq20_sorted.bam"

# --- map (single clean command; backslashes have NO trailing spaces) ---
if [[ "$DATATYPE" == "pe" ]]; then
  bowtie2 -q --phred33 --${MAPPINGPRESET} -p "${SLURM_CPUS_PER_TASK}" -I 0 -X 1500 --fr \
    --rg-id "$SAMPLE_UNIQ_ID" --rg "SM:$SAMPLE_UNIQ_ID" --rg "LB:$SAMPLE_UNIQ_ID" \
    --rg "PU:$PU" --rg "PL:ILLUMINA" \
    -x "$REFBASENAME" \
    -1 "${SAMPLETOMAP}${FASTQSUFFIX1}" -2 "${SAMPLETOMAP}${FASTQSUFFIX2}" \
    -S "$OUTSAM"
else
  bowtie2 -q --phred33 --${MAPPINGPRESET} -p "${SLURM_CPUS_PER_TASK}" \
    --rg-id "$SAMPLE_UNIQ_ID" --rg "SM:$SAMPLE_UNIQ_ID" --rg "LB:$SAMPLE_UNIQ_ID" \
    --rg "PU:$PU" --rg "PL:ILLUMINA" \
    -x "$REFBASENAME" \
    -U "${SAMPLETOMAP}${FASTQSUFFIX1}" \
    -S "$OUTSAM"
fi

# --- SAM -> mapped-only BAM ---
samtools view -bS -F 4 "$OUTSAM" > "$OUTBAM"
rm -f "$OUTSAM"

# --- MAPQ>=20 filter + coordinate sort (single pipe, fixed) ---
samtools view -h -q 20 "$OUTBAM" | samtools sort -@ "${SLURM_CPUS_PER_TASK}" -o "$SORTED" -

if [ -f "$SORTED" ]; then
  echo "SUCCESS: $SAMPLE_UNIQ_ID"; du -h "$SORTED"
  rm -f "$OUTBAM"
else
  echo "ERROR: failed to create $SORTED"; exit 1
fi
echo "Time: $(date)"
