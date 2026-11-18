#!/bin/bash
#SBATCH --job-name=run_qc
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=03:00:00
#SBATCH --output=logs/p07u_run_qc_%A_%a.log

# =============================================================================
# p07u_run_qc.sh - why do the 2024-run samples (1-63) show ~45% more
# heterozygosity than 2023 samples, raw or amplified? One array task per BAM,
# on the first 100 Mb of chr1:
#   1. read-level QC from samtools stats: error rate, base quality, read length,
#      insert size, duplicates
#   2. heterozygous calls at well-covered sites (10-30 reads):
#        het_per_kb   - het calls per 1,000 well-covered sites
#        ab_skewed    - share of hets whose allele balance is outside 0.3-0.7.
#                       DNA from a second animal (contamination) adds many hets
#                       with a skewed balance (e.g. 15:85); true hets sit near 50:50
#        low_af_per_kb - sites per 1,000 where a second allele sits at 5-25% of
#                       the reads (>= 2 reads): the clearest sign of contamination
#        tstv_het     - transition/transversion ratio of the hets. Real variants
#                       give ~2-2.5; random sequencing errors pull it towards 0.5;
#                       DNA damage (C>T) pushes it up
# Submit from /mnt/parscratch/users/bi4og:
#   bash scripts/p07u_run_qc.sh --list          (writes the sample list, prints N)
#   sbatch --array=1-N scripts/p07u_run_qc.sh
#   bash scripts/p07u_run_qc.sh --summary       (after the array finishes)
# =============================================================================

OUT_DIR="physalia/angsd_results/run_qc"
LIST="physalia/sample_lists/bam_list_run_qc.list"
ALL="physalia/angsd_results/sex_det/window_depth/bam_list.txt"     # the 195 BAMs
REFERENCE="/mnt/parscratch/users/bi4og/ncbi_submission/plongirostris.fasta"
REGION="chr1_RagTag_pilon:1-100000000"

group_of() {   # 2024 run = samples 1-63 ("aT" = amplified tissue DNA); "o" = raw 2023 DNA; others = amplified 2023 DNA
    local id=$1 num=${1%%-*}; num=${num%o}
    if [[ "$id" == *o ]]; then echo "2023_raw"
    elif (( num <= 63 )) && [[ "$id" == *aT ]]; then echo "2024_amplified"
    elif (( num <= 63 )); then echo "2024_run"
    else echo "2023_amplified"; fi
}

if [[ "${1:-}" == "--list" ]]; then
    mkdir -p "$(dirname "$LIST")"
    # 2023 raw ("o") and their amplified pairs, the 2024 run's T samples (re-sampled 2023 animals),
    # their 2023 counterparts, and 12 further 2024-run samples
    { grep -E '/[0-9]+o-[A-Za-z0-9]+o_merged' "$ALL"
      grep -E '/[0-9]+-(CAI6|NS1|NS14|NS30|NS32|NS40|NS41|SI18|SI19|SI59|SI7|SI29|SI14|SI1)_merged' "$ALL"
      grep -E '/(60|61|62|63)-SI[0-9a]+T_merged' "$ALL"
      grep -E '/([1-9]|[1-4][0-9]|5[0-8])-[A-Z]+[0-9]+_merged' "$ALL" | awk 'NR%4==1' | head -12
    } | sort -u > "$LIST"
    echo "$(wc -l < "$LIST") BAMs in $LIST:"
    while read -r b; do id=$(basename "$b" | sed 's/_merged.*//'); echo "  $id  $(group_of "$id")"; done < "$LIST"
    echo "Next: sbatch --array=1-$(wc -l < "$LIST") scripts/p07u_run_qc.sh"
    exit 0
fi

if [[ "${1:-}" == "--summary" ]]; then
    F=( "$OUT_DIR"/*.qc.tsv ); [[ -e "${F[0]}" ]] || { echo "no results in $OUT_DIR"; exit 1; }
    head -1 "${F[0]}" > "$OUT_DIR/run_qc_all.tsv"; for f in "${F[@]}"; do tail -n +1 "$f" | sed 1d; done >> "$OUT_DIR/run_qc_all.tsv"
    column -t "$OUT_DIR/run_qc_all.tsv" 2>/dev/null || cat "$OUT_DIR/run_qc_all.tsv"
    echo; echo "Medians by group:"
    Rscript -e 'd <- read.delim(commandArgs(TRUE)[1]); num <- sapply(d, is.numeric)
                m <- aggregate(d[, num], list(group = d$group), median); m$n <- as.vector(table(d$group)[m$group])
                print(m, row.names = FALSE, digits = 4)' "$OUT_DIR/run_qc_all.tsv"
    exit 0
fi

# ---- array task ----
source ~/.bash_profile
conda activate /mnt/parscratch/users/bi4og/conda_envs/owenspopgen
set -eu
for t in samtools bcftools; do command -v $t >/dev/null || { echo "$t not found: conda install -c bioconda $t"; exit 1; }; done
mkdir -p "$OUT_DIR" logs
BAM=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$LIST")
[[ -n "$BAM" ]] || { echo "no BAM on line ${SLURM_ARRAY_TASK_ID} of $LIST"; exit 1; }
ID=$(basename "$BAM" | sed 's/_merged.*//'); GROUP=$(group_of "$ID")
P="$OUT_DIR/$ID"
echo "$ID ($GROUP)  region $REGION"

# 1. read-level QC
samtools stats -@ 2 "$BAM" "$REGION" > "$P.stats"
sn() { grep -P "^SN\t$1" "$P.stats" | cut -f3 | head -1; }
ERR=$(sn "error rate:"); QUAL=$(sn "average quality:"); LEN=$(sn "average length:")
INS=$(sn "insert size average:"); DUP=$(sn "reads duplicated:"); RD=$(sn "reads mapped:")
DUPF=$(awk -v d="$DUP" -v r="$RD" 'BEGIN{printf "%.4f", (r>0)?d/r:0}')

# 2. sites with 10-30 reads (MAPQ>=20, base quality>=20)
WELL=$(samtools depth -@ 2 -q 20 -Q 20 -r "$REGION" "$BAM" | awk '$3>=10 && $3<=30' | wc -l)

# single-sample calls; hets at 10-30 reads with QUAL>=30
bcftools mpileup -f "$REFERENCE" -r "$REGION" -q 20 -Q 20 -a AD,DP -Ou "$BAM" 2>/dev/null \
  | bcftools call -mv -Ou 2>/dev/null \
  | bcftools view -v snps -m2 -M2 -i 'QUAL>=30 && FMT/DP>=10 && FMT/DP<=30 && GT="het"' -Ou \
  | bcftools query -f '%REF\t%ALT\t[%AD]\n' > "$P.hets.tsv"

# sites with 10-30 reads where a second allele is present at a low fraction (5-25%, >= 2 reads).
# Sequencing errors scatter single reads; DNA from another animal gives repeated reads of one allele.
bcftools mpileup -f "$REFERENCE" -r "$REGION" -q 20 -Q 20 -a AD,DP -Ou "$BAM" 2>/dev/null \
  | bcftools query -i 'FMT/DP>=10 && FMT/DP<=30' -f '[%AD]\n' \
  | awk -F, '{t=0; mx=0; for(i=1;i<=NF;i++){t+=$i; if(i>1 && $i>mx) mx=$i} if(t>0 && mx>=2 && mx/t>=0.05 && mx/t<=0.25) n++}
             END{print n+0}' > "$P.lowaf.txt"
LOWAF=$(cat "$P.lowaf.txt")

awk -F'\t' -v id="$ID" -v lowaf="$LOWAF" -v g="$GROUP" -v well="$WELL" -v err="$ERR" -v q="$QUAL" -v len="$LEN" \
    -v ins="$INS" -v dup="$DUPF" '
  function ts(r,a){ return (r a=="AG"||r a=="GA"||r a=="CT"||r a=="TC") }
  { split($3, ad, ","); t=ad[1]+ad[2]; if (t==0) next; n++; ab=ad[2]/t
    if (ab<0.3 || ab>0.7) sk++; if (ts($1,$2)) nts++; else ntv++ }
  END { printf "sample\tgroup\terror_rate\tavg_base_qual\tread_len\tinsert_size\tdup_frac\twell_covered_sites\thets\thet_per_kb\tab_skewed\ttstv_het\tlow_af_per_kb\n"
        printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%d\t%d\t%.4f\t%.3f\t%.2f\t%.4f\n", id, g, err, q, len, ins, dup, well, n,
               (well>0)?1000*n/well:0, (n>0)?sk/n:0, (ntv>0)?nts/ntv:0, (well>0)?1000*lowaf/well:0 }' "$P.hets.tsv" > "$P.qc.tsv"
cat "$P.qc.tsv"
rm -f "$P.stats"
