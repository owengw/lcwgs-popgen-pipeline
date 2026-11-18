# p07v_af_hist.R - plot the allele-fraction histograms from p07v_af_hist.sh, one line per sample,
# one panel per group. Run from /mnt/parscratch/users/bi4og:  Rscript scripts/p07v_af_hist.R
suppressPackageStartupMessages({ library(readr); library(dplyr); library(ggplot2) })
d <- "physalia/angsd_results/run_qc/af_hist"
f <- list.files(d, pattern = "\\.afhist\\.tsv$", full.names = TRUE)
if (!length(f)) stop("no histograms in ", d)
h <- bind_rows(lapply(f, read_tsv, col_types = "cdd", progress = FALSE))
grp <- function(id) {
  num <- as.integer(sub("o$", "", sub("-.*", "", id)))
  ifelse(grepl("o$", id), "2023 raw",
  ifelse(num <= 63 & grepl("aT$", id), "2024 amplified",
  ifelse(num <= 63, "2024 run", "2023 amplified")))
}
tot <- h %>% group_by(sample) %>% summarise(total = sum(sites), .groups = "drop")
h <- h %>% filter(bin_low >= 0) %>% left_join(tot, by = "sample") %>%
  mutate(per_1000 = 1000 * sites / total, af = bin_low + 0.0125, group = grp(sample))
# summary: share of sites with a second allele at 5-30% vs 40-60% of reads
s <- h %>% group_by(group, sample) %>%
  summarise(low_5_30 = sum(per_1000[af >= 0.05 & af < 0.30]),
            mid_40_60 = sum(per_1000[af >= 0.40 & af <= 0.60]),
            high_90_100 = sum(per_1000[af >= 0.90]), .groups = "drop") %>%
  mutate(across(where(is.numeric), ~ round(.x, 2)))
cat("Sites per 1,000 well-covered sites with the second allele at:\n")
print(s %>% arrange(group, sample), n = Inf)
cat("\nGroup medians:\n")
print(s %>% group_by(group) %>% summarise(n = n(), across(c(low_5_30, mid_40_60, high_90_100), median)), n = Inf)
write_tsv(s, file.path(d, "af_summary.tsv"))
p <- ggplot(h %>% filter(af >= 0.05), aes(af, per_1000, group = sample, colour = group)) +
  geom_line(alpha = 0.6) + facet_wrap(~group, ncol = 2) + scale_y_log10() +
  labs(x = "Fraction of reads carrying the second allele (sites with 10-30 reads)",
       y = "Sites per 1,000 (log scale)",
       title = "Clean single animal: peak at 0.5 and little below 0.3") +
  theme_bw(9) + theme(legend.position = "none")
ggsave(file.path(d, "af_histograms.png"), p, width = 9, height = 6, dpi = 150)
cat("\nPlot:", file.path(d, "af_histograms.png"), "\n")
