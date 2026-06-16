#!/bin/bash

#SBATCH --job-name=p02_index_genome
#SBATCH --output=p02_index_genome.log
#SBATCH --mem=16G
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --time=4:00:00

source ~/.bash_profile
conda activate /mnt/parscratch/users/bi4og/conda_envs/owenspopgen

# Genome file path
GENOME=${1:-genome/Plong_genome_flye.fasta}

src=$PWD
REFERENCE=$src/$GENOME
REFBASENAME="${REFERENCE%.*}"

echo "Indexing reference genome: $REFERENCE"

# Create fasta index
echo "Creating fasta index..."
samtools faidx $REFERENCE

# Create sequence dictionary
echo "Creating sequence dictionary..."
picard CreateSequenceDictionary R=$REFERENCE O=$REFBASENAME'.dict'

# Build bowtie2 index
echo "Building bowtie2 index..."
bowtie2-build --threads 8 $REFERENCE $REFBASENAME

echo "Genome indexing complete!"