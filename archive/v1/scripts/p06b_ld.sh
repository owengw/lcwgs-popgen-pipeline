#!/bin/bash
#SBATCH --job-name=p06b_ld
#SBATCH --cpus-per-task=16
#SBATCH --mem=128G
#SBATCH --time=24:00:00
#SBATCH --output=logs/p06b_ld_%A.log

# Load environment BEFORE set -euo pipefail
source ~/.bash_profile
conda activate /mnt/parscratch/users/bi4og/conda_envs/owenspopgen || true

# Fix library path for compute nodes
export LD_LIBRARY_PATH=/mnt/parscratch/users/bi4og/conda_envs/owenspopgen/lib:${LD_LIBRARY_PATH:-}

set -x
set -o pipefail

NGSLD_PATH=/users/bi4og/ngsLD/ngsLD

# Validate ngsLD exists
[[ ! -f "$NGSLD_PATH" ]] && echo "ERROR: ngsLD not found at $NGSLD_PATH" && exit 1

OUT_DIR=${OUT_DIR:-$PWD/physalia/angsd_results}
THREADS=${THREADS:-16}
MAX_KB_DIST=${MAX_KB_DIST:-2000}
MIN_LD=${MIN_LD:-0.5}
SUBSAMPLE_SIZE=${SUBSAMPLE_SIZE:-100000}
REFERENCE=${REFERENCE:-$PWD/genome/Plong_genome_flye.fasta}
BEAGLE="$OUT_DIR/all_samples.beagle.gz"
MAFS="$OUT_DIR/all_samples.mafs.gz"

[[ ! -f "$BEAGLE" ]] && echo "ERROR: Beagle file not found: $BEAGLE" && exit 1
[[ ! -f "$MAFS" ]] && echo "ERROR: MAF file not found: $MAFS" && exit 1

# Count individuals
N_IND=$(zcat "$BEAGLE" | head -1 | awk '{print (NF-3)/3}')

# Subsample sites
(zcat "$BEAGLE" | head -1; \
 zcat "$BEAGLE" | tail -n +2 | \
 awk -v k=$SUBSAMPLE_SIZE 'BEGIN{srand(42)}
    { if (NR <= k) { lines[NR]=$0 } else { r=int(rand()*NR)+1; if(r<=k) lines[r]=$0 } }
    END{for(i=1;i<=k;i++) if(lines[i]!="") print lines[i]}' | \
 sort -t '_' -k1,1 -k2,2n -k3,3n) | gzip > "$OUT_DIR/all_samples.subsampled.beagle.gz"

# Create position file
zcat "$OUT_DIR/all_samples.subsampled.beagle.gz" | tail -n +2 | \
awk '{
    n=split($1,a,"_"); chr=""; for(i=1;i<n;i++) chr=chr (i>1?"_":"") a[i];
    print chr "\t" a[n]
}' | gzip > "$OUT_DIR/all_samples.subsampled.pos.gz"

# Run ngsLD
$NGSLD_PATH \
    --geno "$OUT_DIR/all_samples.subsampled.beagle.gz" \
    --pos "$OUT_DIR/all_samples.subsampled.pos.gz" \
    --probs \
    --n_ind $N_IND \
    --n_sites $(zcat "$OUT_DIR/all_samples.subsampled.pos.gz" | wc -l) \
    --max_kb_dist $MAX_KB_DIST \
    --n_threads $THREADS \
    --out "$OUT_DIR/all_samples.ld" 2>&1 | tee "$OUT_DIR/all_samples.ld.log"

[[ ${PIPESTATUS[0]} -ne 0 ]] && echo "ERROR: ngsLD failed" && exit 1

# Prune SNPs
tail -n +2 "$OUT_DIR/all_samples.ld" > "$OUT_DIR/all_samples.no_header.ld"

python3 - "$OUT_DIR" "$((MAX_KB_DIST*1000))" "$MIN_LD" <<'PYEOF'
import sys
ld_file = sys.argv[1] + "/all_samples.no_header.ld"
max_dist = int(sys.argv[2])
min_r2 = float(sys.argv[3])
output = sys.argv[1] + "/all_samples_unlinked.id"

ld_pairs=[]
with open(ld_file) as f:
    for line in f:
        parts = line.strip().split()
        if len(parts)<4: continue
        pos1,pos2,dist,r2 = parts[0],parts[1],int(parts[2]),float(parts[3])
        if abs(dist)<=max_dist and r2>=min_r2: ld_pairs.append((pos1,pos2,r2))

site_connections={}
for pos1,pos2,r2 in ld_pairs:
    site_connections[pos1]=site_connections.get(pos1,0)+1
    site_connections[pos2]=site_connections.get(pos2,0)+1

removed=set()
for pos1,pos2,r2 in sorted(ld_pairs,key=lambda x:-x[2]):
    if pos1 in removed or pos2 in removed: continue
    if site_connections.get(pos1,0)>=site_connections.get(pos2,0):
        removed.add(pos1)
    else:
        removed.add(pos2)

all_sites=set(site_connections.keys())
unlinked=sorted(all_sites-removed)

with open(output,'w') as out:
    for site in unlinked: out.write(site+'\n')
PYEOF

[[ ! -f "$OUT_DIR/all_samples_unlinked.id" ]] && echo "ERROR: LD pruning failed" && exit 1

# Create LD-pruned SNP list for ANGSD
Rscript --vanilla - "$OUT_DIR" <<'REOF'
library(tidyverse)
args <- commandArgs(trailingOnly=TRUE)
out_dir <- args[1]

pruned_ids <- read_lines(file.path(out_dir,"all_samples_unlinked.id"))
pruned_df <- tibble(id=pruned_ids) %>%
  separate(id,into=c("chromo","position"),sep=":",convert=TRUE)

mafs <- read_tsv(file.path(out_dir,"all_samples.mafs.gz"),
                 col_types=cols(),show_col_types=FALSE) %>%
  select(chromo,position,major,minor)

pruned_sites <- mafs %>% inner_join(pruned_df,by=c("chromo","position"))
write_tsv(pruned_sites,file.path(out_dir,"LDpruned_snps.list"),col_names=FALSE)
REOF

# Generate contig restriction file
cut -f1 "$OUT_DIR/LDpruned_snps.list" | sort -u | \
    awk 'NR==FNR{order[$1]=NR; next} ($1 in order){print order[$1], $1}' \
    "$REFERENCE.fai" - | sort -k1,1n | cut -d' ' -f2 > "$OUT_DIR/LDpruned_contigs.txt"
echo "✓ LDpruned_contigs.txt created ($(wc -l < "$OUT_DIR/LDpruned_contigs.txt") contigs)"

echo "Indexing LD-pruned sites for ANGSD..."
angsd sites index "$OUT_DIR/LDpruned_snps.list"