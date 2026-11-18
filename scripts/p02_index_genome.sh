#!/bin/bash
#SBATCH --job-name=p02_index_genome
#SBATCH --output=logs/p02_index_genome_%j.log
#SBATCH --mem=16G
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --time=4:00:00

# EDIT: default genome now the annotated RagTag/pilon assembly.
# Usage: sbatch scripts/p02_index_genome.sh [genome_fasta]

source ~/.bash_profile
conda activate /mnt/parscratch/users/bi4og/conda_envs/owenspopgen
set -eu   # no pipefail

GENOME=${1:-ncbi_submission/plongirostris.fsa}
src=$PWD
REFERENCE=$src/$GENOME
REFBASENAME="${REFERENCE%.*}"
mkdir -p logs

echo "Indexing reference genome: $REFERENCE"
samtools faidx "$REFERENCE"
picard CreateSequenceDictionary R="$REFERENCE" O="${REFBASENAME}.dict"
bowtie2-build --threads 8 "$REFERENCE" "$REFBASENAME"
echo "Genome indexing complete: $REFBASENAME"
