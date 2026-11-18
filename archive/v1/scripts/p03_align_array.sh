#!/bin/bash 
#SBATCH --job-name=p03_align_array 
#SBATCH --output=p03_align_array_%A_%a.log 
#SBATCH --mem=16G 
#SBATCH --nodes=1 
#SBATCH --ntasks-per-node=1 
#SBATCH --cpus-per-task=8 
#SBATCH --time=24:00:00 
#SBATCH --array=1-333%20 # Adjust to number of lines in your prefix list source 
~/.bash_profile 
conda activate /mnt/parscratch/users/bi4og/conda_envs/owenspopgen 
helpFunction() { 
  echo "" 
  echo "Usage: $0 -g <genome_fasta> -l <prefix_list> -p <adapter_clipped_path>" 
  echo "" 
  echo " -g Genome fasta file (e.g., genome/Plong_genome_flye.fasta)" 
  echo " -l File containing prefixes (one per line, same as used in p01)" 
  echo " -p Path to adapter-clipped fastq files (default: physalia/adapter_clipped)" 
  exit 1 
} 

# Parse options 
while getopts "g:l:p:h" opt; do 
  case "$opt" in 
    g) GENOME="$OPTARG" ;; 
    l) PREFIXLIST="$OPTARG" ;; 
    p) ADAPT_CLIP_PATH="$OPTARG" ;; 
    h) helpFunction ;; 
    *) helpFunction ;; 
  esac 
done 

# Defaults 
GENOME=${GENOME:-genome/Plong_genome_flye.fasta} 
PREFIXLIST=${PREFIXLIST:-physalia/prefix_list.txt} 
ADAPT_CLIP_PATH=${ADAPT_CLIP_PATH:-physalia/adapter_clipped} 
src=$PWD 

mkdir -p "$src/physalia/aligned" 

# Reference setup 

REFERENCE="$src/$GENOME" 
REFBASENAME="${REFERENCE%.*}" 
REFNAME=$(basename "$REFBASENAME") 

# Get prefix for this array task 
PREFIX=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$PREFIXLIST")
if [ -z "$PREFIX" ]; then 
  echo "ERROR: No prefix found for array task $SLURM_ARRAY_TASK_ID" 
  exit 1 
fi 

SAMPLE_UNIQ_ID=$(basename "$PREFIX") 
SAMPLETOMAP="$src/$ADAPT_CLIP_PATH/$SAMPLE_UNIQ_ID" 
SAMPLEBAM="$src/physalia/aligned/$SAMPLE_UNIQ_ID" 

FASTQSUFFIX1=_f_paired.fastq.gz 
FASTQSUFFIX2=_r_paired.fastq.gz 

# Detect paired-end or single-end 
if [[ -f "$SAMPLETOMAP$FASTQSUFFIX2" ]]; then 
  DATATYPE="pe" 
else 
  DATATYPE="se" 
fi 

echo "Processing: $SAMPLE_UNIQ_ID" 
echo "Datatype: $DATATYPE" 
echo "Array ID: $SLURM_ARRAY_TASK_ID" 
echo "Time: $(date)"
 
# Platform unit (just use prefix as identifier) 
PU=$SAMPLE_UNIQ_ID 

# Check input files 
if [[ "$DATATYPE" == "pe" ]]; then 
  for f in "$SAMPLETOMAP$FASTQSUFFIX1" "$SAMPLETOMAP$FASTQSUFFIX2"; do 
    [ -f "$f" ] || { echo "ERROR: File not found: $f"; exit 1; } 
  done 
else 
  [ -f "$SAMPLETOMAP$FASTQSUFFIX1" ] || { echo "ERROR: File not found: $SAMPLETOMAP$FASTQSUFFIX1"; exit 1; } 
fi 

# Map reads 
MAPPINGPRESET=very-sensitive 
if [[ "$DATATYPE" == "pe" ]]; then 
  bowtie2 -q --phred33 --$MAPPINGPRESET -p ${SLURM_CPUS_PER_TASK} -I 0 -X 1500 --fr \ 
  --rg-id $SAMPLE_UNIQ_ID \ 
  --rg SM:$SAMPLE_UNIQ_ID \ 
  --rg LB:$SAMPLE_UNIQ_ID \ 
  --rg PU:$PU \ 
  --rg PL:ILLUMINA \ 
  -x $REFBASENAME \ 
  -1 $SAMPLETOMAP$FASTQSUFFIX1 \ 
  -2 $SAMPLETOMAP$FASTQSUFFIX2 \ 
  -S ${SAMPLEBAM}_${DATATYPE}_bt2_${REFNAME}.sam 
else 
  bowtie2 -q --phred33 --$MAPPINGPRESET -p ${SLURM_CPUS_PER_TASK} \ 
  --rg-id $SAMPLE_UNIQ_ID \ 
  --rg SM:$SAMPLE_UNIQ_ID \ 
  --rg LB:$SAMPLE_UNIQ_ID \ 
  --rg PU:$PU \ 
  --rg PL:ILLUMINA \ 
  -x $REFBASENAME \ 
  -U $SAMPLETOMAP$FASTQSUFFIX1 \ 
  -S ${SAMPLEBAM}_${DATATYPE}_bt2_${REFNAME}.sam 
fi 

# Convert SAM -> BAM (mapped reads only) 
samtools view -bS -F 4 ${SAMPLEBAM}_${DATATYPE}_bt2_${REFNAME}.sam \ 
  > ${SAMPLEBAM}_${DATATYPE}_bt2_${REFNAME}.bam 
  
rm -f ${SAMPLEBAM}_${DATATYPE}_bt2_${REFNAME}.sam 

# Filter & sort (MAPQ =20) 
samtools view -h -q 20 
${SAMPLEBAM}_${DATATYPE}_bt2_${REFNAME}.bam | \ 
  samtools view -buS - | \ samtools sort -@ ${SLURM_CPUS_PER_TASK} -o ${SAMPLEBAM}_${DATATYPE}_bt2_${REFNAME}_minq20_sorted.bam 
  
# Final check 
if [ -f "${SAMPLEBAM}_${DATATYPE}_bt2_${REFNAME}_minq20_sorted.bam" ]; then
  echo "SUCCESS: Completed $SAMPLE_UNIQ_ID" 
  du -h "${SAMPLEBAM}_${DATATYPE}_bt2_${REFNAME}_minq20_sorted.bam" 
else 
  echo "ERROR: Failed to create final sorted BAM file" 
  exit 1 
fi 

echo "Time: $(date)"