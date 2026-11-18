# =============================================================================
# Inbreeding visualisation
# Produces:
#   1. Individual F_HET by population (violin + boxplot + jitter)
#   2. Population mean F_HET with SD (bar chart)
#   3. Combined panel
# All outputs as PDF and SVG
#
# Excludes:
#   - COI, SP, DH (n=3/1, F_HET unreliable at small n)
#   - SI43, SI45 (failed samples, zero heterozygosity)
#   - o-suffix and T-suffix duplicates/originals
# =============================================================================

library(dplyr)
library(tidyr)
library(readr)
library(ggplot2)
library(ggrepel)
library(patchwork)

# =============================================================================
# USER SETTINGS
# =============================================================================

base_dir   <- file.path(Sys.getenv("USERPROFILE"),
                        "OneDrive - MMU/Documents/PhD Analysis/Bioinformatics/Physalia")
fhet_file  <- file.path(base_dir, "inbreeding",
                        "F_HET_individual_corrected.txt")
out_dir    <- file.path(base_dir, "results", "inbreeding")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# Populations shown — COI, SP, DH excluded (n too small for reliable F_HET)
pop_order   <- c("CAI", "NS", "SB", "SI")
pop_colours <- c(
  CAI = "#E69F00",
  NS  = "#56B4E9",
  SB  = "#009E73",
  SI  = "#F0E442"
)

# Samples to exclude
EXCLUDE_SAMPLES  <- c("200-SI43", "216-SI45")
EXCLUDE_POPS     <- c("COI", "SP", "DH")

save_plot <- function(plot, name, width = 8, height = 6) {
  ggsave(file.path(out_dir, paste0(name, ".pdf")), plot,
         width = width, height = height)
  ggsave(file.path(out_dir, paste0(name, ".svg")), plot,
         width = width, height = height, device = "svg")
  cat("Saved:", name, "\n")
}

# =============================================================================
# 1. LOAD F_HET DATA
# =============================================================================

cat("Loading individual F_HET data...\n")

if (!file.exists(fhet_file)) {
  stop("F_HET file not found: ", fhet_file,
       "\nRun p07c_inbreeding.sh first.")
}

fhet_raw <- read_tsv(fhet_file, col_types = cols(), show_col_types = FALSE)

cat("Loaded", nrow(fhet_raw), "individuals\n")

# Extract clean sample ID (e.g. "200-SI43" -> "SI43") and population
fhet <- fhet_raw %>%
  mutate(
    clean_id = str_extract(sample_id, "[A-Z]+[0-9]+[a-z]*$"),
    pop      = str_remove(site, "^pop_map_")
  ) %>%
  # Remove failed samples
  filter(!sample_id %in% EXCLUDE_SAMPLES) %>%
  # Remove duplicates and originals (o-suffix, T-suffix)
  filter(!grepl("o_merged|T_merged", sample_id)) %>%
  # Remove small-n populations
  filter(!pop %in% EXCLUDE_POPS) %>%
  filter(!is.na(pop), pop %in% pop_order) %>%
  mutate(pop = factor(pop, levels = pop_order))

cat("Retained", nrow(fhet), "individuals after filtering\n")
print(count(fhet, pop))

# Flag outliers (> 3 SD from population mean)
fhet <- fhet %>%
  group_by(pop) %>%
  mutate(
    pop_mean_F   = mean(F_HET, na.rm = TRUE),
    pop_sd_F     = sd(F_HET,   na.rm = TRUE),
    is_outlier   = abs(F_HET - pop_mean_F) > 3 * pop_sd_F
  ) %>%
  ungroup()

cat("\nOutlier individuals (>3 SD from pop mean F_HET):\n")
outliers <- filter(fhet, is_outlier)
if (nrow(outliers) > 0) {
  print(select(outliers, sample_id, pop, F_HET))
} else {
  cat("  None\n")
}

# =============================================================================
# 2. POPULATION SUMMARY
# =============================================================================

pop_summary <- fhet %>%
  group_by(pop) %>%
  summarise(
    n              = n(),
    mean_F         = mean(F_HET,   na.rm = TRUE),
    sd_F           = sd(F_HET,     na.rm = TRUE),
    median_F       = median(F_HET, na.rm = TRUE),
    se_F           = sd_F / sqrt(n),
    n_above_0.10   = sum(F_HET > 0.10, na.rm = TRUE),
    n_above_0.25   = sum(F_HET > 0.25, na.rm = TRUE),
    pct_above_0.10 = round(100 * n_above_0.10 / n, 1),
    pct_above_0.25 = round(100 * n_above_0.25 / n, 1),
    .groups        = "drop"
  )

cat("\n=== Population F_HET summary ===\n")
print(pop_summary)

write_tsv(pop_summary,
          file.path(out_dir, "F_HET_population_summary.tsv"))

# =============================================================================
# 3. PLOT — Individual F_HET by population
# =============================================================================

# Reference lines for biological interpretation
f_thresholds <- data.frame(
  yintercept = c(0, 0.125, 0.25),
  label      = c("Panmictic expectation",
                 "1st cousin mating equivalent",
                 "Half-sibling mating equivalent"),
  linetype   = c("solid", "dashed", "dotted")
)

p_individual <- ggplot(fhet,
                       aes(x = pop, y = F_HET,
                           fill = pop, colour = pop)) +
  # Reference lines
  geom_hline(yintercept = 0,     linetype = "solid",
             colour = "grey60",  linewidth = 0.5) +
  geom_hline(yintercept = 0.125, linetype = "dashed",
             colour = "grey60",  linewidth = 0.4) +
  geom_hline(yintercept = 0.25,  linetype = "dotted",
             colour = "grey60",  linewidth = 0.4) +
  # Data
  geom_violin(alpha = 0.25, linewidth = 0.4, trim = FALSE) +
  geom_boxplot(width = 0.12, alpha = 0.8, linewidth = 0.5,
               outlier.shape = NA, colour = "grey30") +
  geom_jitter(width = 0.08, size = 1.2, alpha = 0.45, shape = 16) +
  # Label outliers
  geom_text_repel(
    data        = filter(fhet, is_outlier),
    aes(label   = clean_id),
    size        = 2.5, colour = "grey30",
    box.padding = 0.4, max.overlaps = 20,
    show.legend = FALSE
  ) +
  scale_fill_manual(values = pop_colours,   guide = "none") +
  scale_colour_manual(values = pop_colours, guide = "none") +
  # Annotation for reference lines
  annotate("text", x = 4.6, y = 0.02,    label = "F = 0",
           size = 2.8, colour = "grey50", hjust = 0) +
  annotate("text", x = 4.6, y = 0.145,   label = "F = 0.125",
           size = 2.8, colour = "grey50", hjust = 0) +
  annotate("text", x = 4.6, y = 0.27,    label = "F = 0.25",
           size = 2.8, colour = "grey50", hjust = 0) +
  labs(
    title    = "Individual inbreeding coefficients (F\u2091\u2090\u209c)",
    subtitle = paste0("F = 1 \u2212 (H\u2092\u2095\u209b / H\u2091\u2093\u2090)",
                      "  |  H\u2092\u2095\u209b from 1-sample SFS,",
                      " H\u2091\u2093\u2090 from population \u03c0\n",
                      "Excludes SI43 & SI45 (failed samples)",
                      "  |  COI, SP, DH excluded (n \u2264 3)"),
    x        = "Population",
    y        = expression(F[HET])
  ) +
  coord_cartesian(clip = "off") +
  theme_bw(base_size = 12) +
  theme(
    panel.grid.minor = element_blank(),
    plot.margin      = margin(5, 50, 5, 5)
  )

save_plot(p_individual, "F_HET_individual_by_population",
          width = 9, height = 7)

# =============================================================================
# 4. PLOT — Population mean F_HET
# =============================================================================

p_pop_mean <- ggplot(pop_summary,
                     aes(x = pop, y = mean_F, fill = pop)) +
  geom_hline(yintercept = 0, linetype = "solid",
             colour = "grey60", linewidth = 0.5) +
  geom_col(width = 0.6, colour = "grey40", linewidth = 0.3, alpha = 0.85) +
  geom_errorbar(
    aes(ymin = mean_F - se_F, ymax = mean_F + se_F),
    width = 0.18, linewidth = 0.6, colour = "grey30"
  ) +
  geom_text(
    aes(label = paste0("n=", n, "\n",
                       round(mean_F, 3))),
    vjust  = ifelse(pop_summary$mean_F < 0, 1.5, -0.3),
    size   = 3, colour = "grey20"
  ) +
  scale_fill_manual(values = pop_colours, guide = "none") +
  labs(
    title    = "Population mean inbreeding coefficient",
    subtitle = "Mean F\u2091\u2090\u209c \u00b1 SE  |  COI, SP, DH excluded (n \u2264 3)",
    x        = "Population",
    y        = expression(mean~F[HET])
  ) +
  theme_bw(base_size = 12) +
  theme(panel.grid.minor = element_blank())

save_plot(p_pop_mean, "F_HET_population_mean", width = 7, height = 5)

# =============================================================================
# 5. PLOT — % individuals above inbreeding thresholds
# =============================================================================

threshold_long <- pop_summary %>%
  select(pop, `F > 0.10` = pct_above_0.10,
         `F > 0.25` = pct_above_0.25) %>%
  pivot_longer(cols      = starts_with("F"),
               names_to  = "threshold",
               values_to = "pct") %>%
  mutate(threshold = factor(threshold,
                            levels = c("F > 0.10", "F > 0.25")))

p_thresholds <- ggplot(threshold_long,
                       aes(x = pop, y = pct,
                           fill = pop, alpha = threshold)) +
  geom_col(position = position_dodge(width = 0.7),
           width = 0.65, colour = "grey40", linewidth = 0.3) +
  geom_text(
    aes(label = paste0(round(pct, 1), "%")),
    position = position_dodge(width = 0.7),
    vjust    = -0.4, size = 2.8, colour = "grey20"
  ) +
  scale_fill_manual(values = pop_colours, guide = "none") +
  scale_alpha_manual(
    values = c("F > 0.10" = 0.85, "F > 0.25" = 0.45),
    name   = "Threshold"
  ) +
  labs(
    title    = "Proportion of individuals above inbreeding thresholds",
    subtitle = paste0("Darker = F > 0.10 (moderate)  |",
                      "  Lighter = F > 0.25 (high, equivalent to half-sibling mating)"),
    x        = "Population",
    y        = "Individuals (%)"
  ) +
  theme_bw(base_size = 12) +
  theme(
    panel.grid.minor = element_blank(),
    legend.position  = "bottom"
  )

save_plot(p_thresholds, "F_HET_threshold_proportions", width = 8, height = 5)

# =============================================================================
# 6. COMBINED INBREEDING PANEL
# =============================================================================

p_combined <- (p_individual / (p_pop_mean | p_thresholds)) +
  plot_layout(heights = c(1.6, 1)) +
  plot_annotation(
    title    = "Population inbreeding analysis",
    subtitle = paste0("F\u2091\u2090\u209c = 1 \u2212 (H\u2092\u2095\u209b / H\u2091\u2093\u2090)  |  ",
                      "H\u2091\u2093\u2090 estimated from population \u03c0 (corrected thetaStat)  |  ",
                      "Populations CAI, NS, SB, SI only"),
    theme    = theme(
      plot.title    = element_text(size = 14, face = "bold"),
      plot.subtitle = element_text(size = 10)
    )
  )

save_plot(p_combined, "inbreeding_combined_panel", width = 14, height = 14)

cat("\nAll inbreeding plots saved to:", out_dir, "\n")

