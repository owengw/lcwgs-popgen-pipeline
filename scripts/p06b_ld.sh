#!/bin/bash
#SBATCH --job-name=p06b_ld
#SBATCH --cpus-per-task=16
#SBATCH --mem=128G
#SBATCH --time=24:00:00
#SBATCH --output=logs/p06b_ld_%A.log

# Load environment BEFORE set -eu
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
REFERENCE=${REFERENCE:?ERROR: export REFERENCE (assembly the BAMs were mapped to)}
BEAGLE="$OUT_DIR/all_samples.beagle.gz"
MAFS="$OUT_DIR/all_samples.mafs.gz"

[[ ! -f "$BEAGLE" ]] && echo "ERROR: Beagle file not found: $BEAGLE" && exit 1
[[ ! -f "$MAFS" ]] && echo "ERROR: MAF file not found: $MAFS" && exit 1

# Count individuals
N_IND=$(zcat "$BEAGLE" | head -1 | awk '{print (NF-3)/3}')

# Subsample sites
# FIX: the old sort used `-t '_' -k1,1 -k2,2n -k3,3n`, but markers look like
# chr1_RagTag_pilon_100012300 (pos is the LAST underscore field, not field 2/3),
# so rows came out in LEXICOGRAPHIC order and ngsLD died with
# "invalid distance between adjacent sites". We now build an explicit
# chr <TAB> pos <TAB> line key, sort numerically on position, then strip the key.
# This is also robust to unplaced contigs (contig_123_RagTag_pilon_456).
(zcat "$BEAGLE" | head -1; \
 zcat "$BEAGLE" | tail -n +2 | \
 awk -v k=$SUBSAMPLE_SIZE 'BEGIN{srand(42)}
    { if (NR <= k) { lines[NR]=$0 } else { r=int(rand()*NR)+1; if(r<=k) lines[r]=$0 } }
    END{for(i=1;i<=k;i++) if(lines[i]!="") print lines[i]}' | \
 awk 'BEGIN{OFS="\t"}{ n=split($1,a,"_"); pos=a[n]; chr=a[1];
        for(i=2;i<n;i++) chr=chr "_" a[i];
        print chr, pos, $0 }' | \
 sort -k1,1 -k2,2n | cut -f3-) | gzip > "$OUT_DIR/all_samples.subsampled.beagle.gz"

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
# FIXES:
#  - system python3 / Rscript+tidyverse are unusable here (GLIBC_2.28 missing,
#    R not on PATH). Rewritten in awk: no external dependencies.
#  - a stale Flye REFERENCE default silently produced a Flye-contig SNP list.
tail -n +2 "$OUT_DIR/all_samples.ld" > "$OUT_DIR/all_samples.no_header.ld"

MAXD=$((MAX_KB_DIST*1000))

awk -v maxd="$MAXD" -v minr2="$MIN_LD" '
{
    if (NF < 4) next
    p1=$1; p2=$2; d=$3+0; r2=$4+0
    if (d<0) d=-d
    if (d<=maxd && r2>=minr2) { n++; A[n]=p1; B[n]=p2; R[n]=r2; conn[p1]++; conn[p2]++ }
}
END{
    for(i=1;i<=n;i++) ord[i]=i
    for(i=2;i<=n;i++){ v=ord[i]; j=i-1
        while(j>=1 && R[ord[j]]<R[v]){ ord[j+1]=ord[j]; j-- }
        ord[j+1]=v }
    for(k=1;k<=n;k++){ i=ord[k]
        if(A[i] in removed || B[i] in removed) continue
        if(conn[A[i]]>=conn[B[i]]) removed[A[i]]=1; else removed[B[i]]=1 }
    for(s in conn) if(!(s in removed)) print s
}' "$OUT_DIR/all_samples.no_header.ld" | sort > "$OUT_DIR/all_samples_unlinked.id"

[[ -s "$OUT_DIR/all_samples_unlinked.id" ]] || { echo "ERROR: LD pruning produced no sites"; exit 1; }
echo "Unlinked sites retained: $(wc -l < "$OUT_DIR/all_samples_unlinked.id")"

# ANGSD sites file: chromo <TAB> pos <TAB> major <TAB> minor
zcat "$MAFS" | tail -n +2 | awk -v ids="$OUT_DIR/all_samples_unlinked.id" '
BEGIN{ while((getline line < ids) > 0) keep[line]=1 }
{ if(($1":"$2) in keep) print $1"\t"$2"\t"$3"\t"$4 }' > "$OUT_DIR/LDpruned_snps.list"

[[ -s "$OUT_DIR/LDpruned_snps.list" ]] || { echo "ERROR: LDpruned_snps.list empty"; exit 1; }
echo "LD-pruned SNPs: $(wc -l < "$OUT_DIR/LDpruned_snps.list")"

[[ -f "$REFERENCE.fai" ]] || { echo "ERROR: .fai not found: $REFERENCE.fai"; exit 1; }
cut -f1 "$OUT_DIR/LDpruned_snps.list" | sort -u | \
    awk 'NR==FNR{order[$1]=NR; next} ($1 in order){print order[$1], $1}' \
    "$REFERENCE.fai" - | sort -k1,1n | cut -d' ' -f2 > "$OUT_DIR/LDpruned_contigs.txt"
echo "LDpruned_contigs.txt created ($(wc -l < "$OUT_DIR/LDpruned_contigs.txt") contigs)"

echo "Indexing LD-pruned sites for ANGSD..."
angsd sites index "$OUT_DIR/LDpruned_snps.list"