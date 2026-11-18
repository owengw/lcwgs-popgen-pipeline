# =============================================================================
# ROH visualisation — ROHan output
# Plots individual and population-level ROH metrics from ROHan
#
# Key metrics:
#   theta     — heterozygosity outside ROH (genome-wide diversity)
#   pct_roh   — % genome in ROH (F_ROH proxy)
#   avg_roh_kb — average ROH length in kb
#   total_roh_bp — total bp in ROH
#
# IMPORTANT CAVEAT:
#   ROH detection is limited by the fragmented assembly (~4289 contigs).
#   ROH spanning contig boundaries are missed, so F_ROH values are
#   conservative underestimates. Most individuals showing zero ROH likely
#   have genuine ROH that are split across contigs.
#
# Populations: CAI, NS, SB, SI
# =============================================================================

library(dplyr)
library(tidyr)
library(readr)
library(ggplot2)
library(patchwork)

# =============================================================================
# USER SETTINGS
# =============================================================================

base_dir <- file.path(Sys.getenv("USERPROFILE"),
                      "OneDrive - MMU/Documents/PhD Analysis/Bioinformatics/Physalia")
roh_dir  <- file.path(base_dir, "roh")
out_dir  <- file.path(base_dir, "results", "roh")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

pop_order   <- c("CAI", "NS", "SB", "SI")
pop_colours <- c(CAI = "#E69F00", NS = "#56B4E9",
                 SB  = "#009E73", SI = "#F0E442")

save_plot <- function(plot, name, width = 10, height = 7) {
  ggsave(file.path(out_dir, paste0(name, ".pdf")), plot,
         width = width, height = height)
  ggsave(file.path(out_dir, paste0(name, ".svg")), plot,
         width = width, height = height, device = "svg")
  cat("Saved:", name, "\n")
}

# =============================================================================
# 1. LOAD DATA
# =============================================================================

ind_roh <- read.table(file.path(roh_dir, "individual_ROH_all.tsv"),
                      header = TRUE, sep = "\t", stringsAsFactors = FALSE) %>%
  filter(site %in% pop_order) %>%
  mutate(
    site          = factor(site, levels = pop_order),
    total_roh_Mb  = n_roh / 1e6,
    has_roh       = pct_roh > 0
  )

pop_roh <- read.table(file.path(roh_dir, "population_ROH_summary.tsv"),
                      header = TRUE, sep = "\t", stringsAsFactors = FALSE) %>%
  filter(population %in% pop_order) %>%
  mutate(population = factor(population, levels = pop_order))

cat("Individuals loaded:", nrow(ind_roh), "\n")
cat("Individuals with detectable ROH:", sum(ind_roh$has_roh), "\n")
print(count(ind_roh, site, has_roh))

# =============================================================================
# 2. PLOT 1 — % genome in ROH per individual (dot plot)
# =============================================================================

# Most individuals have zero — show as strip with jitter
# Highlight individuals with ROH > 0

p_pct_roh <- ggplot(ind_roh,
                    aes(x = site, y = pct_roh,
                        colour = site, fill = site)) +
  geom_hline(yintercept = 0, colour = "grey80", linewidth = 0.4) +
  geom_jitter(data   = filter(ind_roh, !has_roh),
              width  = 0.25, size = 1.5, alpha = 0.3, shape = 16) +
  geom_jitter(data   = filter(ind_roh, has_roh),
              width  = 0.15, size = 3, alpha = 0.9, shape = 18) +
  geom_boxplot(width = 0.15, alpha = 0, linewidth = 0.6,
               outlier.shape = NA, colour = "grey30") +
  scale_colour_manual(values = pop_colours, guide = "none") +
  scale_fill_manual(values   = pop_colours, guide = "none") +
  labs(
    title    = "Proportion of genome in ROH per individual",
    subtitle = paste0(
      "Diamonds = individuals with detectable ROH  |  ",
      "Circles = zero ROH detected\n",
      "\u26a0 Conservative estimate: ROH spanning contig boundaries undetectable ",
      "(~4,289 contigs)"
    ),
    x = "Population",
    y = "Genome in ROH (%)"
  ) +
  theme_bw(base_size = 12) +
  theme(
    panel.grid.minor = element_blank(),
    plot.subtitle    = element_text(size = 9, colour = "grey40")
  )

save_plot(p_pct_roh, "roh_pct_genome", width = 8, height = 6)

# =============================================================================
# 3. PLOT 2 — Individual theta (heterozygosity outside ROH)
# =============================================================================

p_theta <- ggplot(ind_roh,
                  aes(x = site, y = theta,
                      colour = site, fill = site)) +
  geom_violin(alpha = 0.2, linewidth = 0.4, trim = TRUE) +
  geom_jitter(width = 0.15, size = 1.2, alpha = 0.5, shape = 16) +
  geom_boxplot(width = 0.1, alpha = 0.8, linewidth = 0.5,
               outlier.shape = NA, colour = "grey30") +
  scale_colour_manual(values = pop_colours, guide = "none") +
  scale_fill_manual(values   = pop_colours, guide = "none") +
  labs(
    title    = "Individual heterozygosity outside ROH (ROHan \u03b8)",
    subtitle = paste0("Genome-wide heterozygosity estimated outside ROH windows  |  ",
                      "Comparable to individual het from 1-sample SFS"),
    x = "Population",
    y = expression(theta~"(outside ROH)")
  ) +
  theme_bw(base_size = 12) +
  theme(
    panel.grid.minor = element_blank(),
    plot.subtitle    = element_text(size = 9, colour = "grey40")
  )

save_plot(p_theta, "roh_theta", width = 8, height = 6)

# =============================================================================
# 4. PLOT 3 — Population summary comparison
# =============================================================================

pop_long <- pop_roh %>%
  select(population, mean_pct_roh, mean_avg_roh_kb) %>%
  pivot_longer(cols      = c(mean_pct_roh, mean_avg_roh_kb),
               names_to  = "metric",
               values_to = "value") %>%
  mutate(metric = recode(metric,
                         mean_pct_roh    = "Mean % genome in ROH",
                         mean_avg_roh_kb = "Mean avg. ROH length (kb)"))

p_pop_summary <- ggplot(pop_long,
                        aes(x = population, y = value,
                            fill = population)) +
  geom_col(width = 0.6, colour = "grey30", linewidth = 0.3) +
  geom_text(aes(label = round(value, 3)),
            vjust = -0.4, size = 3.2, colour = "grey20") +
  scale_fill_manual(values = pop_colours, guide = "none") +
  facet_wrap(~metric, scales = "free_y") +
  labs(
    title    = "Population-level ROH summary",
    subtitle = paste0("Mean values per population  |  ",
                      "Conservative estimates due to fragmented assembly"),
    x = "Population",
    y = "Value"
  ) +
  theme_bw(base_size = 12) +
  theme(
    panel.grid.minor = element_blank(),
    strip.text       = element_text(face = "bold"),
    plot.subtitle    = element_text(size = 9, colour = "grey40")
  )

save_plot(p_pop_summary, "roh_population_summary", width = 10, height = 6)

# =============================================================================
# 5. PLOT 4 — ROH vs F_HET comparison
# Load F_HET data and compare with ROH theta
# =============================================================================

fhet_file <- file.path(base_dir, "inbreeding",
                       "F_HET_individual_corrected.txt")

if (file.exists(fhet_file)) {
  fhet <- read.table(fhet_file, header = TRUE, sep = "\t",
                     stringsAsFactors = FALSE) %>%
    filter(population %in% pop_order) %>%
    mutate(population = factor(population, levels = pop_order))

  # Match on sample ID
  combined <- ind_roh %>%
    inner_join(fhet, by = c("sample_id", "site" = "population"))

  if (nrow(combined) > 0) {
    p_compare <- ggplot(combined,
                        aes(x = F_HET, y = pct_roh,
                            colour = site, shape = has_roh)) +
      geom_point(size = 2.5, alpha = 0.7) +
      scale_colour_manual(values = pop_colours, name = "Population") +
      scale_shape_manual(values = c("FALSE" = 16, "TRUE" = 18),
                         labels = c("FALSE" = "No ROH detected",
                                    "TRUE"  = "ROH detected"),
                         name   = NULL) +
      labs(
        title    = "ROH (%genome) vs F_HET per individual",
        subtitle = paste0("Both measure inbreeding but via different approaches  |  ",
                          "Low concordance expected due to assembly fragmentation"),
        x = expression(F[HET]~"(1-sample SFS method)"),
        y = "% genome in ROH (ROHan)"
      ) +
      theme_bw(base_size = 12) +
      theme(
        panel.grid.minor = element_blank(),
        plot.subtitle    = element_text(size = 9, colour = "grey40"),
        legend.position  = "bottom"
      )

    save_plot(p_compare, "roh_vs_fhet", width = 9, height = 7)
  } else {
    cat("NOTE: Could not match ROH and F_HET individuals — skipping comparison plot\n")
  }
} else {
  cat("NOTE: F_HET file not found — skipping comparison plot\n")
  cat("Expected:", fhet_file, "\n")
}

# =============================================================================
# 6. COMBINED PANEL
# =============================================================================

p_combined <- (p_theta | p_pct_roh) /
              p_pop_summary +
  plot_annotation(
    title    = "Runs of Homozygosity analysis (ROHan)",
    subtitle = paste0(
      "Individual ROH estimated from lcWGS using hidden Markov model  |  ",
      "CAI, NS, SB, SI  |  n=193 individuals\n",
      "\u26a0 Results are conservative: fragmented assembly (~4,289 contigs) ",
      "prevents detection of ROH spanning contig boundaries"
    ),
    theme = theme(
      plot.title    = element_text(size = 14, face = "bold"),
      plot.subtitle = element_text(size = 9, colour = "grey40")
    )
  ) +
  plot_layout(heights = c(1, 1))

save_plot(p_combined, "roh_combined_panel", width = 14, height = 12)

# =============================================================================
# 7. SUMMARY TABLE
# =============================================================================

cat("\n=== ROH summary ===\n")
cat("Individuals with detectable ROH by population:\n")
print(
  ind_roh %>%
    group_by(site) %>%
    summarise(
      n_total        = n(),
      n_with_roh     = sum(has_roh),
      pct_with_roh   = round(100 * mean(has_roh), 1),
      mean_pct_roh   = round(mean(pct_roh), 4),
      mean_theta     = round(mean(theta, na.rm = TRUE), 6),
      .groups = "drop"
    )
)

cat("\n=== Interpretation notes ===\n")
cat("1. Most individuals show zero detectable ROH — largely an assembly artefact\n")
cat("2. ROHan theta (outside ROH) is comparable to individual het from 1-sample SFS\n")
cat("3. Population ranking by mean_pct_roh should match F_HET ranking if both\n")
cat("   are detecting the same inbreeding signal\n")
cat("4. Discordance between ROH and F_HET may reflect:\n")
cat("   - Assembly fragmentation limiting ROH detection\n")
cat("   - F_HET captures ancient inbreeding; ROH captures recent inbreeding\n")
cat("   - Different genomic scales: F_HET genome-wide, ROH requires long tracts\n")
cat("\nAll plots saved to:", out_dir, "\n")
