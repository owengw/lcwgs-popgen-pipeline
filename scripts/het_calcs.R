
library(dplyr)
library(tidyr)
library(readr)
library(ggplot2)
library(patchwork)

base_dir   <- "C:/Users/23608589/OneDrive - MMU/Documents/PhD Analysis/Bioinformatics/Physalia"
het_dir    <- file.path(base_dir, "heterozygosity")
fst_dir    <- file.path(base_dir, "fst")
pca_dir    <- file.path(base_dir, "PCA")
meta_file  <- file.path(base_dir, "metadata.tsv")
out_dir    <- file.path(base_dir, "results")
dir.create(out_dir, showWarnings = FALSE)

populations <- c("CAI", "NS", "SB", "SI")
pop_colours <- c(CAI = "#E69F00", NS = "#56B4E9",
                 SB = "#009E73", SI = "#F0E442")

save_plot <- function(plot, name, width = 10, height = 6) {
  ggsave(file.path(out_dir, paste0(name, ".pdf")), plot,
         width = width, height = height)
  ggsave(file.path(out_dir, paste0(name, ".svg")), plot,
         width = width, height = height, device = "svg")
  cat("Saved:", name, "\n")
}

# ---------------------------------------------------------------------------
# Load individual heterozygosity summary
# ---------------------------------------------------------------------------
ind_het <- read_tsv(
  file.path(het_dir, "het_all_individuals_individual_summary.txt"),
  show_col_types = FALSE
) %>%
  mutate(
    sample_id = sub("_merged_pe_bt2.*", "", sample_id),
    Ho_ind    = mean_heterozygosity   # posterior heterozygosity per individual
  )

# Remove duplicates
ind_het <- ind_het %>%
  filter(!grepl("(aT|[^a]T|o)$", sample_id))

# Join population info from metadata
metadata <- read_tsv(meta_file, show_col_types = FALSE)

ind_het <- ind_het %>%
  left_join(metadata %>% select(sample_id, population), by = "sample_id") %>%
  filter(population %in% populations) %>%
  mutate(population = factor(population, levels = populations))

cat("Individuals after filtering:", nrow(ind_het), "\n")
cat("Per population:\n")
print(count(ind_het, population))

# ---------------------------------------------------------------------------
# Calculate He per population from mafs (allele frequencies)
# He = mean of 2pq across all SNP sites
# ---------------------------------------------------------------------------
read_gz <- function(path, ...) {
  read_tsv(gzcon(file(path, "rb")), show_col_types = FALSE, ...)
}

he_pop <- lapply(populations, function(pop) {
  f <- file.path(fst_dir, paste0("pop_map_", pop, ".mafs.gz"))
  if (!file.exists(f)) { warning("mafs not found: ", pop); return(NULL) }
  mafs <- read_gz(f) %>%
    mutate(He_site = 2 * knownEM * (1 - knownEM))
  tibble(
    population = pop,
    He         = mean(mafs$He_site, na.rm = TRUE),
    He_sd      = sd(mafs$He_site,   na.rm = TRUE),
    n_sites_He = nrow(mafs)
  )
}) %>%
  bind_rows() %>%
  mutate(population = factor(population, levels = populations))

cat("\nHe per population (from mafs):\n")
print(he_pop)

# ---------------------------------------------------------------------------
# Ho summary per population from individual values
# ---------------------------------------------------------------------------
ho_pop <- ind_het %>%
  group_by(population) %>%
  summarise(
    n_ind  = n(),
    Ho     = mean(Ho_ind, na.rm = TRUE),
    Ho_sd  = sd(Ho_ind,   na.rm = TRUE),
    .groups = "drop"
  )

# Combine
ho_he_summary <- ho_pop %>%
  left_join(he_pop, by = "population") %>%
  mutate(Fis = 1 - (Ho / He))

cat("\nFinal Ho/He summary:\n")
print(ho_he_summary %>% mutate(across(where(is.numeric), ~round(., 5))))

write_tsv(ho_he_summary, file.path(out_dir, "Ho_He_summary.tsv"))
write_tsv(ind_het %>% select(sample_id, population, Ho_ind, n_sites,
                             observed_het_proportion),
          file.path(out_dir, "Ho_He_individual.tsv"))

# ---------------------------------------------------------------------------
# Plots
# ---------------------------------------------------------------------------

# Grouped bar: Ho vs He
ho_he_long <- ho_he_summary %>%
  select(population, Ho, He, Ho_sd, He_sd) %>%
  pivot_longer(c(Ho, He), names_to = "metric", values_to = "value") %>%
  mutate(sd = if_else(metric == "Ho", Ho_sd, He_sd))

p_bar <- ggplot(ho_he_long,
                aes(x = population, y = value, fill = metric)) +
  geom_col(position = position_dodge(0.7), width = 0.6, alpha = 0.85) +
  geom_errorbar(aes(ymin = value - sd, ymax = value + sd),
                position = position_dodge(0.7),
                width = 0.2, linewidth = 0.4) +
  scale_fill_manual(values = c(Ho = "#E69F00", He = "#0072B2"),
                    labels = c(Ho = "Observed (Ho)",
                               He = "Expected (He)"),
                    name = NULL) +
  labs(title    = "Observed vs expected heterozygosity per population",
       subtitle = paste0("Ho = mean posterior individual heterozygosity ",
                         "(~181k sites)\n",
                         "He = mean 2pq from SNP allele frequencies"),
       x = NULL, y = "Heterozygosity") +
  theme_bw(base_size = 12) +
  theme(legend.position = "top")

p_fis <- ggplot(ho_he_summary,
                aes(x = population, y = Fis, fill = population)) +
  geom_col(width = 0.6, alpha = 0.85) +
  geom_hline(yintercept = 0, linetype = "dashed",
             colour = "grey40", linewidth = 0.5) +
  scale_fill_manual(values = pop_colours) +
  labs(title    = "Fixation index (Fis) per population",
       subtitle = "Fis = 1 - (Ho/He)",
       x = NULL, y = "Fis") +
  theme_bw(base_size = 12) +
  theme(legend.position = "none")

p_violin <- ggplot(ind_het,
                   aes(x = population, y = Ho_ind,
                       fill = population, colour = population)) +
  geom_violin(alpha = 0.35, linewidth = 0.4, trim = FALSE) +
  geom_boxplot(width = 0.12, alpha = 0.8, outlier.shape = NA,
               linewidth = 0.4, colour = "grey30") +
  geom_point(data = ho_he_summary,
             aes(x = population, y = He),
             shape = 23, size = 3, fill = "white",
             colour = "grey20", stroke = 1,
             inherit.aes = FALSE) +
  scale_fill_manual(values   = pop_colours) +
  scale_colour_manual(values = pop_colours) +
  scale_x_discrete(labels = function(x) {
    n <- ho_pop$n_ind[match(x, as.character(ho_pop$population))]
    paste0(x, "\n(n=", n, ")")
  }) +
  labs(title    = "Individual observed heterozygosity (Ho) per population",
       subtitle = "Distribution of per-individual Ho; diamond = population He",
       x = NULL, y = "Observed heterozygosity (Ho)") +
  theme_bw(base_size = 12) +
  theme(legend.position      = "none",
        panel.grid.major.x   = element_blank())

save_plot(p_bar,    "Ho_He_barplot",        width = 8,  height = 5)
save_plot(p_fis,    "Fis_barplot",          width = 7,  height = 5)
save_plot(p_violin, "Ho_individual_violin", width = 8,  height = 5)

p_combined <- (p_bar | p_fis) / p_violin +
  plot_annotation(
    title = "Heterozygosity analysis",
    theme = theme(plot.title = element_text(size = 14, face = "bold"))
  )
save_plot(p_combined, "heterozygosity_combined", width = 14, height = 10)
cat("Heterozygosity analysis complete.\n")

# Read heterozygosity summary and extract population-level rows
het_summary <- read_tsv(
  file.path(het_dir, "heterozygosity_summary.txt"),
  show_col_types = FALSE
) %>%
  filter(population %in% c("pop_map_CAI", "pop_map_NS",
                           "pop_map_SB",  "pop_map_SI")) %>%
  mutate(
    population = sub("pop_map_", "", population),
    population = factor(population, levels = populations),
    # theta_pi / n_segregating_sites gives per-site diversity
    pi_per_site = theta_pi / n_segregating_sites
  ) %>%
  rename(
    Ho            = heterozygosity,
    n_seg_sites   = n_segregating_sites
  ) %>%
  select(population, n_individuals, n_seg_sites, Ho, pi_per_site)

# Join individual-level Ho for violin plot
ind_het_clean <- read_tsv(
  file.path(het_dir, "het_all_individuals_individual_summary.txt"),
  show_col_types = FALSE
) %>%
  mutate(sample_id = sub("_merged_pe_bt2.*", "", sample_id)) %>%
  filter(!grepl("(aT|[^a]T|o)$", sample_id)) %>%
  left_join(read_tsv(meta_file, show_col_types = FALSE) %>%
              select(sample_id, population),
            by = "sample_id") %>%
  filter(population %in% populations) %>%
  mutate(population = factor(population, levels = populations))

print(het_summary)
