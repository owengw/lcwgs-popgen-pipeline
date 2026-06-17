# =============================================================================
# HWE visualisation
# Produces:
#   1. % significant sites per population (bar chart, three thresholds)
#   2. Per-site F distribution by population (violin + boxplot)
#   3. Combined panel
#   4. Supplementary LD decay figure with caveat
# All outputs as PDF and SVG
#
# Populations: CAI, NS, SB, SI
# Data from: hwe/pop_map_{POP}_hwe_results.tsv.gz
#            hwe/hwe_summary_all_populations.tsv
# =============================================================================

library(dplyr)
library(tidyr)
library(readr)
library(ggplot2)
library(patchwork)

# =============================================================================
# USER SETTINGS
# =============================================================================

base_dir  <- file.path(Sys.getenv("USERPROFILE"),
                        "OneDrive - MMU/Documents/PhD Analysis/Bioinformatics/Physalia")
hwe_dir   <- file.path(base_dir, "hwe")
ld_dir    <- file.path(base_dir, "ld decay")
out_dir   <- file.path(base_dir, "results", "hwe")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

pop_order   <- c("CAI", "NS", "SB", "SI")
pop_colours <- c(
  CAI = "#E69F00",
  NS  = "#56B4E9",
  SB  = "#009E73",
  SI  = "#F0E442"
)

save_plot <- function(plot, name, width = 8, height = 6) {
  ggsave(file.path(out_dir, paste0(name, ".pdf")), plot,
         width = width, height = height)
  ggsave(file.path(out_dir, paste0(name, ".svg")), plot,
         width = width, height = height, device = "svg")
  cat("Saved:", name, "\n")
}

# =============================================================================
# 1. LOAD HWE SUMMARY
# =============================================================================

cat("Loading HWE summary...\n")

summary_file <- file.path(hwe_dir, "hwe_summary_all_populations.tsv")
if (!file.exists(summary_file)) {
  stop("HWE summary not found: ", summary_file,
       "\nRun p08c_hwe.sh first.")
}

hwe_summary <- read_tsv(summary_file, col_types = cols(),
                         show_col_types = FALSE) %>%
  filter(population %in% pop_order) %>%
  mutate(population = factor(population, levels = pop_order))

cat("Populations loaded:", paste(hwe_summary$population, collapse = ", "), "\n")
print(hwe_summary)

# =============================================================================
# 2. LOAD PER-SITE HWE RESULTS
# =============================================================================

cat("\nLoading per-site HWE results...\n")

persite_list <- lapply(pop_order, function(pop) {
  f <- file.path(hwe_dir, paste0("pop_map_", pop, "_hwe_results.tsv.gz"))
  if (!file.exists(f)) {
    cat("WARNING: Per-site file not found:", f, "\n")
    return(NULL)
  }
  dat <- read_tsv(f, col_types = cols(), show_col_types = FALSE)
  dat$population <- pop
  dat
})

persite <- bind_rows(persite_list[!sapply(persite_list, is.null)]) %>%
  mutate(population = factor(population, levels = pop_order))

cat("Total sites loaded:", nrow(persite), "\n")
print(count(persite, population))

# =============================================================================
# 3. PLOT — % significant sites by population
# =============================================================================

# Reshape for grouped bar chart showing three thresholds
sig_long <- hwe_summary %>%
  select(population,
         `p < 0.05`           = pct_sig_p005,
         `p < 0.001`          = n_sig_p001,
         `Bonferroni`         = n_sig_bonferroni) %>%
  mutate(
    `p < 0.001`  = round(100 * `p < 0.001`  / hwe_summary$n_sites_tested, 2),
    `Bonferroni` = round(100 * `Bonferroni` / hwe_summary$n_sites_tested, 2)
  ) %>%
  pivot_longer(cols      = c(`p < 0.05`, `p < 0.001`, `Bonferroni`),
               names_to  = "threshold",
               values_to = "pct") %>%
  mutate(threshold = factor(threshold,
                             levels = c("p < 0.05", "p < 0.001", "Bonferroni")))

# Expected % under true HWE at p<0.05 = 5% (false positive rate)
p_sig <- ggplot(sig_long,
                aes(x = population, y = pct,
                    fill = population, alpha = threshold)) +
  geom_hline(yintercept = 5, linetype = "dashed",
             colour = "grey50", linewidth = 0.5) +
  annotate("text", x = -Inf, y = 5.8,
           label = "5% expected under HWE",
           hjust = 0, size = 3, colour = "grey50") +
  geom_col(position = position_dodge(width = 0.75),
           width = 0.7, colour = "grey40", linewidth = 0.3) +
  geom_text(
    aes(label = paste0(round(pct, 1), "%")),
    position = position_dodge(width = 0.75),
    vjust = -0.4, size = 2.8, colour = "grey20"
  ) +
  scale_fill_manual(values = pop_colours, guide = "none") +
  scale_alpha_manual(
    values = c("p < 0.05" = 0.90, "p < 0.001" = 0.60, "Bonferroni" = 0.35),
    name   = "Significance threshold"
  ) +
  scale_y_continuous(limits = c(0, max(sig_long$pct) * 1.15),
                     labels = function(x) paste0(x, "%")) +
  labs(
    title    = "HWE deviation per population",
    subtitle = paste0("% of SNPs significantly deviating from Hardy-Weinberg equilibrium\n",
                      "Dashed line = 5% expected under true HWE (false positive rate)"),
    x        = "Population",
    y        = "SNPs deviating from HWE (%)"
  ) +
  theme_bw(base_size = 12) +
  theme(
    panel.grid.minor = element_blank(),
    legend.position  = "bottom",
    legend.key.size  = unit(0.4, "cm")
  )

save_plot(p_sig, "hwe_pct_significant", width = 8, height = 6)

# =============================================================================
# 4. PLOT — Per-site F distribution by population
# =============================================================================

# Flag outliers > 3SD
persite <- persite %>%
  group_by(population) %>%
  mutate(
    pop_mean_F = mean(F_site,   na.rm = TRUE),
    pop_sd_F   = sd(F_site,     na.rm = TRUE),
    is_outlier = abs(F_site - pop_mean_F) > 3 * pop_sd_F
  ) %>%
  ungroup()

p_fsite <- ggplot(persite,
                  aes(x = population, y = F_site,
                      fill = population, colour = population)) +
  geom_hline(yintercept = 0, linetype = "solid",
             colour = "grey60", linewidth = 0.5) +
  geom_violin(alpha = 0.25, linewidth = 0.4, trim = TRUE) +
  geom_boxplot(width = 0.12, alpha = 0.8, linewidth = 0.5,
               outlier.shape = NA, colour = "grey30") +
  scale_fill_manual(values = pop_colours,   guide = "none") +
  scale_colour_manual(values = pop_colours, guide = "none") +
  labs(
    title    = "Per-site inbreeding coefficient from HWE analysis",
    subtitle = paste0("F\u209b\u1d35\u209c\u1d49 = 1 \u2212 (H\u2092\u2095\u209b / H\u2091\u2093\u2090) per SNP  |  ",
                      "Positive values indicate heterozygote deficit"),
    x        = "Population",
    y        = expression(F[site])
  ) +
  theme_bw(base_size = 12) +
  theme(panel.grid.minor = element_blank())

save_plot(p_fsite, "hwe_fsite_distribution", width = 8, height = 6)

# =============================================================================
# 5. COMBINED HWE PANEL
# =============================================================================

p_hwe_combined <- (p_sig | p_fsite) +
  plot_annotation(
    title    = "Hardy-Weinberg Equilibrium analysis",
    subtitle = paste0("GL-based HWE test (HardyWeinberg R package)  |  ",
                      "SNPs at p < 1\u00d710\u207b\u2076  |  CAI, NS, SB, SI only"),
    theme    = theme(
      plot.title    = element_text(size = 14, face = "bold"),
      plot.subtitle = element_text(size = 10)
    )
  )

save_plot(p_hwe_combined, "hwe_combined_panel", width = 14, height = 6)

# =============================================================================
# 6. SUPPLEMENTARY — LD decay with caveat
# =============================================================================

ld_file <- file.path(ld_dir, "ld_decay_all_populations.tsv")

if (!file.exists(ld_file)) {
  cat("\nNOTE: LD decay file not found — skipping supplementary LD figure\n")
} else {
  cat("\nGenerating supplementary LD decay figure...\n")

  ld_decay <- read_tsv(ld_file, col_types = cols(), show_col_types = FALSE) %>%
    filter(population %in% pop_order) %>%
    mutate(population = factor(population, levels = pop_order))

  # Only plot bins with sufficient pairs
  ld_decay_filt <- ld_decay %>%
    filter(n_pairs >= 3)

  p_ld <- ggplot(ld_decay_filt,
                 aes(x = dist_kb, y = mean_r2,
                     colour = population, group = population)) +
    geom_line(linewidth = 0.8, alpha = 0.8) +
    geom_point(aes(size = n_pairs), alpha = 0.6, shape = 16) +
    scale_colour_manual(values = pop_colours, name = "Population") +
    scale_size_continuous(name   = "SNP pairs\nper bin",
                          range  = c(1, 4),
                          breaks = c(5, 10, 20, 50)) +
    scale_x_continuous(labels = function(x) paste0(x, "kb")) +
    labs(
      title    = "LD decay — supplementary",
      subtitle = paste0(
        "Mean r\u00b2 between SNP pairs binned by physical distance (1kb windows)\n",
        "\u26a0 CAUTION: Fragmented assembly (~4,289 contigs) limits LD decay interpretation.\n",
        "Only within-contig pairs computed; few pairs per distance bin produces\n",
        "noisy estimates. Results should be interpreted with caution."
      ),
      x = "Distance (kb)",
      y = expression(Mean~r^2)
    ) +
    theme_bw(base_size = 12) +
    theme(
      panel.grid.minor  = element_blank(),
      plot.subtitle     = element_text(size = 9, colour = "grey40"),
      legend.position   = "right"
    )

  supp_dir <- file.path(base_dir, "results", "supplementary")
  dir.create(supp_dir, showWarnings = FALSE, recursive = TRUE)

  ggsave(file.path(supp_dir, "supplementary_ld_decay.pdf"), p_ld,
         width = 10, height = 6)
  ggsave(file.path(supp_dir, "supplementary_ld_decay.svg"), p_ld,
         width = 10, height = 6, device = "svg")
  cat("Saved: supplementary_ld_decay\n")
}

cat("\nAll HWE plots saved to:", out_dir, "\n")

