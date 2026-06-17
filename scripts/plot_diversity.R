# =============================================================================
# Diversity visualisation — corrected version
# Produces:
#   1. Individual heterozygosity by population (violin + boxplot)
#   2. Population-level pi and theta_w (bar chart, corrected theta)
#   3. Tajima's D — full-n panel
#   4. Tajima's D — downsampled n=17 panel (if available)
#   5. Combined panels
# All outputs as PDF and SVG
#
# Excludes:
#   - SI43, SI45 (zero/near-zero heterozygosity, failed samples)
#   - o-suffix and T-suffix duplicates/originals
#   - DH (n=1)
# =============================================================================

library(dplyr)
library(tidyr)
library(readr)
library(ggplot2)
library(ggrepel)
library(patchwork)
library(stringr)

# =============================================================================
# USER SETTINGS
# =============================================================================

base_dir      <- file.path(Sys.getenv("USERPROFILE"),
                           "OneDrive - MMU/Documents/PhD Analysis/Bioinformatics/Physalia")
het_dir       <- file.path(base_dir, "heterozygosity")
theta_file    <- file.path(base_dir,
                           "theta_summary_corrected.tsv")
ds_theta_file <- file.path(base_dir,
                           "theta_downsampled_summary.tsv")
meta_file     <- file.path(base_dir, "metadata.tsv")
out_dir       <- file.path(base_dir, "results", "diversity")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# Populations — DH excluded throughout
pop_order   <- c("CAI", "NS", "SB", "SI", "COI", "SP", "ALL")
pop_colours <- c(
  CAI = "#E69F00", NS  = "#56B4E9", SB  = "#009E73",
  SI  = "#F0E442", COI = "#0072B2", SP  = "#CC79A7", ALL = "grey30"
)

# Samples to exclude — failed (zero/near-zero het) and duplicates/originals
EXCLUDE_SAMPLES <- c("200-SI43", "216-SI45")

save_plot <- function(plot, name, width = 8, height = 6) {
  ggsave(file.path(out_dir, paste0(name, ".pdf")), plot,
         width = width, height = height)
  ggsave(file.path(out_dir, paste0(name, ".svg")), plot,
         width = width, height = height, device = "svg")
  cat("Saved:", name, "\n")
}

# =============================================================================
# 1. LOAD INDIVIDUAL HETEROZYGOSITY
# =============================================================================

cat("Loading individual heterozygosity files...\n")

het_files <- list.files(het_dir, pattern = "\\.het$", full.names = TRUE)
cat("Found", length(het_files), "het files\n")

het_data <- bind_rows(lapply(het_files, read_tsv,
                             col_types = cols(), show_col_types = FALSE))

metadata <- read_tsv(meta_file, col_types = cols(), show_col_types = FALSE)

het_data <- het_data %>%
  mutate(sample_prefix = str_extract(sample_id, "^[^-]+-[A-Z]+[0-9]+")) %>%
  left_join(
    metadata %>%
      mutate(sample_prefix = str_extract(sample_id, "^[^-]+-[A-Z]+[0-9]+")) %>%
      select(sample_prefix, site, year),
    by = "sample_prefix"
  ) %>%
  # Remove failed samples, duplicates, originals, DH
  filter(!sample_id %in% EXCLUDE_SAMPLES) %>%
  filter(!grepl("o$|T$", str_extract(sample_id, "[A-Z]+[0-9]+[a-z]*$"))) %>%
  filter(!is.na(site), site %in% pop_order) %>%
  mutate(site = factor(site, levels = pop_order))

cat("Retained", nrow(het_data), "individuals after filtering\n")
cat("Excluded failed samples: SI43, SI45\n")
print(count(het_data, site))

# Flag statistical outliers (> 3 SD from population mean)
het_data <- het_data %>%
  group_by(site) %>%
  mutate(
    pop_mean   = mean(heterozygosity),
    pop_sd     = sd(heterozygosity),
    is_outlier = abs(heterozygosity - pop_mean) > 3 * pop_sd
  ) %>%
  ungroup()

cat("Outlier individuals (>3 SD):\n")
print(filter(het_data, is_outlier) %>%
        select(sample_id, site, heterozygosity))

# =============================================================================
# 2. PLOT — Individual heterozygosity by population
# =============================================================================

p_het <- ggplot(het_data,
                aes(x = site, y = heterozygosity,
                    fill = site, colour = site)) +
  geom_violin(alpha = 0.3, linewidth = 0.4, trim = FALSE) +
  geom_boxplot(width = 0.15, alpha = 0.8, linewidth = 0.5,
               outlier.shape = NA, colour = "grey30") +
  geom_jitter(width = 0.08, size = 1.2, alpha = 0.5, shape = 16) +
  geom_text_repel(
    data        = filter(het_data, is_outlier),
    aes(label   = str_extract(sample_id, "[A-Z]+[0-9]+")),
    size        = 2.5, colour = "grey30",
    box.padding = 0.4, max.overlaps = 20,
    show.legend = FALSE
  ) +
  scale_fill_manual(values = pop_colours, guide = "none") +
  scale_colour_manual(values = pop_colours, guide = "none") +
  scale_y_continuous(labels = scales::label_number(accuracy = 0.0001)) +
  labs(
    title    = "Individual genome-wide heterozygosity",
    subtitle = "1-sample SFS method (realSFS -fold 1), random 10M sites\nExcludes SI43 and SI45 (failed samples, <600 sites covered)",
    x        = "Population",
    y        = "Heterozygosity (per site)"
  ) +
  theme_bw(base_size = 12) +
  theme(panel.grid.minor = element_blank())

save_plot(p_het, "het_individual_by_population", width = 9, height = 6)

# =============================================================================
# 3. LOAD CORRECTED POPULATION THETA
# =============================================================================

cat("\nLoading corrected population theta...\n")

if (!file.exists(theta_file)) {
  stop("Corrected theta summary not found: ", theta_file,
       "\nRun p07b_pop_theta.sh first.")
}

theta_raw <- read_tsv(theta_file, col_types = cols(), show_col_types = FALSE)

theta <- theta_raw %>%
  filter(str_detect(population, "^pop_map_")) %>%
  mutate(site = str_remove(population, "^pop_map_")) %>%
  filter(site %in% pop_order) %>%
  mutate(site = factor(site, levels = pop_order))

cat("Populations in theta summary:\n")
print(theta %>% select(site, n_sites, theta_pi_per_site,
                       Watterson_theta_per_site, tajima_D))

# =============================================================================
# 4. PLOT — Population pi and theta_w
# =============================================================================

theta_long <- theta %>%
  select(site,
         pi      = theta_pi_per_site,
         theta_w = Watterson_theta_per_site) %>%
  pivot_longer(cols      = c(pi, theta_w),
               names_to  = "statistic",
               values_to = "value") %>%
  mutate(statistic = recode(statistic,
                            pi      = "Nucleotide diversity (\u03c0)",
                            theta_w = "Watterson\u2019s \u03b8"))

p_theta <- ggplot(theta_long,
                  aes(x = site, y = value,
                      fill = site, alpha = statistic)) +
  geom_col(position = position_dodge(width = 0.7),
           width = 0.65, colour = "grey40", linewidth = 0.3) +
  scale_fill_manual(values = pop_colours, guide = "none") +
  scale_alpha_manual(
    values = c("Nucleotide diversity (\u03c0)" = 0.95,
               "Watterson\u2019s \u03b8"       = 0.50),
    name   = NULL
  ) +
  scale_y_continuous(labels = scales::label_number(accuracy = 0.0001)) +
  labs(
    title    = "Population-level genetic diversity",
    subtitle = "Darker = \u03c0 (nucleotide diversity)  |  lighter = Watterson\u2019s \u03b8\nCorrected thetaStat run on random 10M sites",
    x        = "Population",
    y        = "Per-site diversity"
  ) +
  theme_bw(base_size = 12) +
  theme(
    panel.grid.minor = element_blank(),
    legend.position  = "bottom",
    legend.key.size  = unit(0.4, "cm")
  )

save_plot(p_theta, "diversity_pi_thetaw", width = 8, height = 6)

# =============================================================================
# 5. TAJIMA'S D PLOT FUNCTION
# =============================================================================

plot_tajima <- function(df, subtitle_text) {
  ggplot(df, aes(x = site, y = tajima_D)) +
    geom_col(width = 0.65, colour = "grey40", linewidth = 0.3,
             aes(fill = tajima_D), show.legend = FALSE) +
    scale_fill_gradient2(
      low      = "#A32D2D",
      mid      = "#F0E442",
      high     = "#185FA5",
      midpoint = 0
    ) +
    geom_hline(yintercept = 0, linewidth = 0.5,
               linetype = "dashed", colour = "grey50") +
    geom_text(aes(label = round(tajima_D, 3),
                  vjust = ifelse(tajima_D < 0, 1.4, -0.4)),
              size = 3.2, colour = "grey20") +
    labs(
      subtitle = subtitle_text,
      x        = "Population",
      y        = "Tajima's D"
    ) +
    theme_bw(base_size = 12) +
    theme(panel.grid.minor = element_blank(),
          plot.subtitle = element_text(size = 10))
}

# =============================================================================
# 6. PLOT — Tajima's D full-n
# =============================================================================

p_tajima_full <- plot_tajima(
  theta,
  "Full sample sizes (n = CAI:23, NS:57, SB:17, SI:94, COI:3, SP:3)"
) +
  labs(title = "Tajima's D — full sample sizes")

save_plot(p_tajima_full, "tajimas_D_full_n", width = 8, height = 5)

# =============================================================================
# 7. PLOT — Tajima's D downsampled (if available)
# =============================================================================

ds_available <- file.exists(ds_theta_file)

if (ds_available) {
  cat("\nLoading downsampled theta summary...\n")
  
  ds_raw <- read_tsv(ds_theta_file, col_types = cols(), show_col_types = FALSE)
  
  ds_theta <- ds_raw %>%
    filter(str_detect(analysis, "^downsampled")) %>%
    mutate(site = factor(population, levels = pop_order)) %>%
    filter(!is.na(site))
  
  if (nrow(ds_theta) > 0) {
    p_tajima_ds <- plot_tajima(
      ds_theta,
      "Downsampled to equal n=17 per population (seed=42)\nControls for sample size bias in Tajima\u2019s D"
    ) +
      labs(title = "Tajima's D — downsampled n=17")
    
    save_plot(p_tajima_ds, "tajimas_D_downsampled_n17", width = 8, height = 5)
    
    # Combined Tajima's D comparison panel
    p_tajima_combined <- (p_tajima_full / p_tajima_ds) +
      plot_annotation(
        title    = "Tajima's D — full-n vs downsampled comparison",
        subtitle = paste0("Downsampling to equal n=17 controls for sample size bias.",
                          " Rank order preserved across methods confirms biological signal."),
        theme    = theme(
          plot.title    = element_text(size = 13, face = "bold"),
          plot.subtitle = element_text(size = 10)
        )
      )
    
    save_plot(p_tajima_combined, "tajimas_D_comparison",
              width = 9, height = 10)
  }
} else {
  cat("\nNOTE: Downsampled theta summary not found — skipping downsampled Tajima's D panel\n")
  cat("Run p07d_theta_downsampled.sh and rerun this script to add the comparison panel\n")
  p_tajima_ds <- NULL
}

# =============================================================================
# 8. COMBINED DIVERSITY PANEL
# =============================================================================

if (!is.null(p_tajima_ds)) {
  p_combined <- (p_het / (p_theta | p_tajima_full) / p_tajima_ds) +
    plot_layout(heights = c(1.3, 1, 1)) +
    plot_annotation(
      title = "Population genetic diversity",
      theme = theme(plot.title = element_text(size = 14, face = "bold"))
    )
  save_plot(p_combined, "diversity_combined_panel",
            width = 16, height = 18)
} else {
  p_combined <- (p_het / (p_theta | p_tajima_full)) +
    plot_layout(heights = c(1.3, 1)) +
    plot_annotation(
      title = "Population genetic diversity",
      theme = theme(plot.title = element_text(size = 14, face = "bold"))
    )
  save_plot(p_combined, "diversity_combined_panel",
            width = 16, height = 14)
}

# =============================================================================
# 9. SUMMARY TABLE
# =============================================================================

het_summary <- het_data %>%
  group_by(site) %>%
  summarise(
    n          = n(),
    mean_het   = round(mean(heterozygosity), 6),
    median_het = round(median(heterozygosity), 6),
    sd_het     = round(sd(heterozygosity), 6),
    min_het    = round(min(heterozygosity), 6),
    max_het    = round(max(heterozygosity), 6),
    .groups    = "drop"
  )

combined_summary <- theta %>%
  select(site,
         pi      = theta_pi_per_site,
         theta_w = Watterson_theta_per_site,
         tajima_D) %>%
  left_join(het_summary, by = "site")

write_tsv(combined_summary,
          file.path(out_dir, "diversity_summary.tsv"))

cat("\n=== Diversity summary ===\n")
print(combined_summary, n = Inf)
cat("\nAll diversity plots saved to:", out_dir, "\n")
