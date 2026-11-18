# p07t_compare.R - do two sequencings of the same animal agree?
#   raw vs amplified (2023 run):  the "o" sample (raw DNA) vs its 2023 amplified original
#   2024 run vs 2023 run:         the T / aT sample (tissue DNA, 2024 run) vs the 2023 sample of that animal
# Run from /mnt/parscratch/users/bi4og:  Rscript scripts/p07t_compare.R
d <- "physalia/angsd_results/het_batch_test"
f <- list.files(d, pattern = "\\.het$", recursive = TRUE, full.names = TRUE)
if (!length(f)) stop("no .het files under ", d)
h <- do.call(rbind, lapply(f, read.delim))
h$animal <- sub("(o|a?T)$", "", sub("^[^-]+-", "", h$sample_id))
h$kind <- ifelse(grepl("o$", h$sample_id), "raw_2023",
          ifelse(grepl("a?T$", h$sample_id), "run_2024", "amplified_2023"))
for (m in sort(unique(h$min_depth))) {
  x <- h[h$min_depth == m, ]
  w <- reshape(x[, c("animal", "kind", "heterozygosity")], idvar = "animal", timevar = "kind", direction = "wide")
  names(w) <- sub("heterozygosity\\.", "", names(w))
  cat(sprintf("\n=== setMinDepthInd %d ===\n", m))
  for (cmp in list(c("raw_2023", "raw (o) vs amplified, both 2023 run"),
                   c("run_2024", "2024 run (T/aT) vs 2023 run, same animal"))) {
    if (!cmp[1] %in% names(w) || !"amplified_2023" %in% names(w)) next
    p <- w[!is.na(w[[cmp[1]]]) & !is.na(w$amplified_2023), c("animal", "amplified_2023", cmp[1])]
    if (!nrow(p)) next
    p$ratio <- round(p[[cmp[1]]] / p$amplified_2023, 3)
    cat(sprintf("\n-- %s: %d animals --\n", cmp[2], nrow(p)))
    print(p, row.names = FALSE)
    if (nrow(p) >= 3) {
      tt <- t.test(log(p[[cmp[1]]]), log(p$amplified_2023), paired = TRUE)
      cat(sprintf("Median ratio %.3f (1 = no difference); paired t-test on log scale p = %.3g\n", median(p$ratio), tt$p.value))
    } else cat("Median ratio", median(p$ratio), "(too few animals for a test)\n")
  }
}
