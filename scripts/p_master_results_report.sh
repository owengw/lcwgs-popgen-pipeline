#!/bin/bash
#SBATCH --job-name=master_report
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=00:15:00
#SBATCH --output=logs/master_report_%j.log

# =============================================================================
# p_master_results_report.sh — assembles ALL corrected results from this
# session into one PDF: whitelisted files only (no directory scanning, so
# nothing stale/oversized/wrong can sneak in — see earlier report-building
# lessons this session). Run p_final_results_plots.sh first.
# =============================================================================

set -x
OUT_DIR=${OUT_DIR:-physalia/angsd_results}
REPORT=physalia/REPORT_FINAL_MASTER
FINAL="$REPORT/Plestiodon_longirostris_full_results.pdf"

rm -rf "$REPORT"
mkdir -p "$REPORT/pages"

# ---- explicit whitelist: source_file | display_title -----------------------
declare -a ITEMS=(
    # --- Population diversity & neutrality ---
    "$OUT_DIR/diversity_plots/01_heterozygosity_by_population.pdf|Genome-wide heterozygosity by population"
    "$OUT_DIR/diversity_plots/05_tajimasD_by_population.pdf|Tajima's D by population (sample-size-corrected)"
    "$OUT_DIR/diversity_plots/03_pcangsd_F_by_population.pdf|Individual inbreeding coefficient (PCAngsd F) by population"
    # --- Population differentiation & structure ---
    "$OUT_DIR/final_plots/01_fst_heatmap.pdf|Pairwise FST between populations"
    "$OUT_DIR/pca_admixture_plots/PCA_PC1_PC2.pdf|PCA: PC1 vs PC2"
    "$OUT_DIR/pca_admixture_plots/PCA_PC1_PC3.pdf|PCA: PC1 vs PC3"
    "$OUT_DIR/pca_admixture_plots/PCAngsd_admix_K2_barplot.pdf|PCAngsd admixture, K=2"
    "$OUT_DIR/pca_admixture_plots/PCAngsd_admix_K3_barplot.pdf|PCAngsd admixture, K=3"
    "$OUT_DIR/pca_admixture_plots/PCAngsd_admix_K4_barplot.pdf|PCAngsd admixture, K=4"
    # --- Inbreeding, ROH, HWE ---
    "$OUT_DIR/final_plots/02_roh_theta_by_population.pdf|ROHan genome-wide theta (outside ROH) by population"
    "$OUT_DIR/final_plots/03_hwe_deviation_by_population.pdf|Hardy-Weinberg deviation by population"
    # --- Relatedness ---
    "$OUT_DIR/final_plots/04_relatedness_by_population.pdf|Relatedness (lcMLkin) by population"
    # --- Sex determination ---
    "$OUT_DIR/sex_det/assignment/plots/sex_ratio_by_site.pdf|Sex ratio by site"
    # --- Sample QC / amplification bias ---
    "$OUT_DIR/final_plots/05_amp_bias_breadth.pdf|Amplification bias: coverage breadth"
    "$OUT_DIR/final_plots/06_amp_bias_gc.pdf|Amplification bias: GC content"
)

echo "=================================================================="
echo "Checking whitelisted files"
echo "=================================================================="
MISSING=0
PRESENT=()
for item in "${ITEMS[@]}"; do
    f="${item%%|*}"; title="${item##*|}"
    if [[ -s "$f" ]]; then
        echo "  OK   : $title"
        PRESENT+=("$item")
    else
        echo "  MISS : $title  ($f)"
        MISSING=$((MISSING+1))
    fi
done
echo ""
echo "Present: ${#PRESENT[@]} / ${#ITEMS[@]}  (missing: $MISSING)"
[[ ${#PRESENT[@]} -eq 0 ]] && { echo "ERROR: nothing to include"; exit 1; }

# ---- cover page ----
module load R/4.4.1-foss-2022b 2>/dev/null || true
Rscript --vanilla -e "
pdf('$REPORT/pages/00_cover.pdf', width=10, height=7)
plot.new()
text(0.5, 0.65, 'Plestiodon longirostris\nfull population genomics results', cex=1.7, font=2)
text(0.5, 0.5, paste('Generated:', Sys.time()), cex=0.85)
text(0.5, 0.4, 'RagTag chromosome-scale assembly; BAQ/mismatch-corrected ANGSD pipeline', cex=0.75, col='grey30')
dev.off()
"

# ---- key summary numbers as a single text page ----
Rscript --vanilla - "$OUT_DIR" "$REPORT" << 'REOF'
args <- commandArgs(trailingOnly = TRUE)
out_dir <- args[1]; report_dir <- args[2]

lines <- c("KEY RESULTS SUMMARY", "")

f <- file.path(out_dir, "theta_downsampled", "theta_downsampled_summary.tsv")
if (file.exists(f)) {
    d <- read.delim(f); d <- d[d$analysis == "downsampled_n16", ]
    lines <- c(lines, "Diversity (n=16, sample-size-corrected):",
        sprintf("  %-6s pi/site=%.5f  Tajima's D=%.3f", d$population, d$theta_pi_per_site, d$tajima_D), "")
}
f <- file.path(out_dir, "fst_corrected")
if (dir.exists(f)) {
    lines <- c(lines, "Pairwise FST:")
    for (ff in list.files(f, pattern="fst.global.txt$", full.names=TRUE)) {
        v <- strsplit(readLines(ff), "\t")[[1]]
        lines <- c(lines, sprintf("  %-10s %s", gsub("\\..*","",basename(ff)), v[2]))
    }
    lines <- c(lines, "")
}
f <- file.path(out_dir, "sex_det", "assignment", "sex_assignments_by_site.tsv")
if (file.exists(f)) {
    d <- read.delim(f)
    lines <- c(lines, "Sex assignment (chr5-based, n per site):",
        sprintf("  %-6s n=%d  male=%d  female=%d  (%.1f%% male)", d$site, d$n_total, d$n_male_het, d$n_female_het, d$pct_male), "")
}

pdf(file.path(report_dir, "pages", "00b_summary.pdf"), width=10, height=8)
plot.new()
text(0.02, 0.95, paste(lines, collapse="\n"), adj=c(0,1), family="mono", cex=0.75)
dev.off()
REOF

# ---- assemble ----
N=1
for item in "${PRESENT[@]}"; do
    f="${item%%|*}"
    cp "$f" "$REPORT/pages/$(printf "%02d" $N)_$(basename "$f")"
    N=$((N+1))
done

command -v gs >/dev/null 2>&1 || { echo "ERROR: ghostscript not found"; exit 1; }
ALL_PAGES=("$REPORT/pages/00_cover.pdf" "$REPORT/pages/00b_summary.pdf")
while IFS= read -r f; do ALL_PAGES+=("$f"); done < <(find "$REPORT/pages" -name "*.pdf" ! -name "00_cover.pdf" ! -name "00b_summary.pdf" | sort)

gs -q -dNOPAUSE -dBATCH -sDEVICE=pdfwrite -sOutputFile="$FINAL" "${ALL_PAGES[@]}"

echo ""
echo "=================================================================="
echo "DONE: $FINAL"
ls -lh "$FINAL"
echo "=================================================================="
