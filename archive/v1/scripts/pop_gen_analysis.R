# =============================================================================
# Population statistics summary table
# Metrics:
#   Ho      — per-site heterozygosity from beagle GLs (heterozygosity_summary)
#   theta_pi — nucleotide diversity from SFS/pestPG (independent of Ho)
#   theta_W  — Watterson's theta from pestPG
#   Tajima D — from pestPG
#   Fst      — pairwise weighted Fst from global fst files
#   Allelic richness — rarefied from SFS files
#   Private alleles  — from SFS files
#   F_HET    — inbreeding coefficient from F_HET_individual.txt
# =============================================================================

library(dplyr)
library(tidyr)
library(readr)
library(ggplot2)
library(knitr)
library(kableExtra)
library(patchwork)

# =============================================================================
# USER SETTINGS
# =============================================================================

base_dir     <- "C:/Users/23608589/OneDrive - MMU/Documents/PhD Analysis/Bioinformatics/Physalia"
het_dir      <- file.path(base_dir, "heterozygosity")
fst_dir      <- file.path(base_dir, "fst")
inbreed_dir  <- file.path(base_dir, "inbreeding")
out_dir      <- file.path(base_dir, "results")
dir.create(out_dir, showWarnings = FALSE)

# All populations for table (include small ones with caveats)
all_populations  <- c("CAI", "NS", "SB", "SI", "COI", "DH", "SP")
# Populations for plots
plot_populations <- c("CAI", "NS", "SB", "SI")

pop_colours <- c(
  CAI = "#E69F00", NS = "#56B4E9", SB = "#009E73", SI = "#F0E442",
  COI = "#0072B2", DH = "#D55E00", SP = "#CC79A7"
)

save_plot <- function(plot, name, width = 10, height = 6) {
  ggsave(file.path(out_dir, paste0(name, ".pdf")), plot,
         width = width, height = height)
  ggsave(file.path(out_dir, paste0(name, ".svg")), plot,
         width = width, height = height, device = "svg")
  cat("Saved:", name, "\n")
}

read_gz <- function(path, ...) {
  read_tsv(gzcon(file(path, "rb")), show_col_types = FALSE, ...)
}

# =============================================================================
# 1. OBSERVED HETEROZYGOSITY (Ho) from heterozygosity_summary.txt
# Ho = mean per-site posterior heterozygosity from beagle GLs
# Calculated across ~181k sites (segregating + invariant)
# =============================================================================

het_raw <- read_tsv(
  file.path(het_dir, "heterozygosity_summary.txt"),
  show_col_types = FALSE
)

# Extract full-population rows only (pop_map_XXX format)
ho_table <- het_raw %>%
  filter(grepl("^pop_map_[A-Z]+$", population)) %>%
  mutate(population = sub("pop_map_", "", population)) %>%
  filter(population %in% all_populations) %>%
  rename(
    n_ind_Ho    = n_individuals,
    n_seg_sites = n_segregating_sites,
    Ho          = heterozygosity
  ) %>%
  select(population, n_ind_Ho, n_seg_sites, Ho) %>%
  mutate(population = factor(population, levels = all_populations))

cat("Ho per population:\n")
print(ho_table)

# =============================================================================
# 2. NUCLEOTIDE DIVERSITY (θπ) AND TAJIMA'S D from pestPG files
# θπ and θW are per-site values (tP/nSites and tW/nSites)
# These are calculated from the SFS independently of Ho
# =============================================================================

theta_list <- lapply(all_populations, function(pop) {
  f <- file.path(fst_dir,
                 paste0("pop_map_", pop, ".thetas.idx.pestPG"))
  if (!file.exists(f)) {
    warning("pestPG not found for ", pop); return(NULL)
  }
  df <- read_tsv(f, show_col_types = FALSE)
  colnames(df)[1] <- "index"
  
  tibble(
    population   = pop,
    n_sites_theta = sum(df$nSites, na.rm = TRUE),
    theta_pi     = sum(df$tP, na.rm = TRUE) /
      sum(df$nSites, na.rm = TRUE),
    theta_W      = sum(df$tW, na.rm = TRUE) /
      sum(df$nSites, na.rm = TRUE),
    Tajima_D     = mean(df$Tajima, na.rm = TRUE)
  )
}) %>%
  bind_rows() %>%
  mutate(population = factor(population, levels = all_populations))

cat("\nTheta per population:\n")
print(theta_list %>% mutate(across(where(is.numeric), ~round(., 6))))

# =============================================================================
# 3. ALLELIC RICHNESS AND PRIVATE ALLELES from SFS files
# Rarefied to smallest population size
# =============================================================================

read_sfs <- function(pop) {
  f <- file.path(fst_dir, paste0("pop_map_", pop, ".sfs"))
  if (!file.exists(f)) return(NULL)
  vals <- as.numeric(strsplit(trimws(readLines(f)), "\\s+")[[1]])
  vals
}

sfs_list <- lapply(all_populations, function(pop) {
  sfs <- read_sfs(pop)
  if (is.null(sfs)) return(NULL)
  n_alleles <- length(sfs) - 1
  n_inds    <- n_alleles / 2
  tibble(
    population      = pop,
    n_inds_sfs      = n_inds,
    seg_sites       = sum(sfs[2:n_alleles]),
    private_alleles = sfs[2],
    sfs             = list(sfs)
  )
}) %>%
  bind_rows()

# Rarefaction for allelic richness
rarefaction_AR <- function(sfs, n_ref) {
  n_total <- length(sfs) - 1
  if (n_total < n_ref) return(NA_real_)
  ar <- sum(sapply(2:length(sfs), function(j) {
    freq   <- round(sfs[j])
    absent <- n_total - freq
    if (freq == 0) return(0)
    absent <- max(absent, 0)
    1 - exp(lchoose(absent, n_ref) - lchoose(n_total, n_ref))
  }))
  ar
}

n_ref <- 2 * min(sfs_list$n_inds_sfs)
cat("\nRarefaction reference:", n_ref, "gene copies\n")

sfs_stats <- sfs_list %>%
  mutate(
    allelic_richness = mapply(
      function(sfs, n) rarefaction_AR(sfs, n_ref), sfs, n_inds_sfs
    ),
    population = factor(population, levels = all_populations)
  ) %>%
  select(population, n_inds_sfs, seg_sites,
         private_alleles, allelic_richness)

# =============================================================================
# 4. INBREEDING (F_HET) from F_HET_individual.txt
# =============================================================================

fhet <- read_tsv(
  file.path(inbreed_dir, "F_HET_individual.txt"),
  show_col_types = FALSE
) %>%
  mutate(
    sample_id  = sub("_merged_pe_bt2.*", "", sample_id),
    population = sub("pop_map_", "", population)
  ) %>%
  filter(
    !grepl("(aT|[^a]T|o)$", sample_id),
    population %in% all_populations
  )

fhet_summary <- fhet %>%
  group_by(population) %>%
  summarise(
    n_ind_fhet = n(),
    mean_F_HET = mean(F_HET, na.rm = TRUE),
    sd_F_HET   = sd(F_HET,   na.rm = TRUE),
    n_inbred   = sum(F_HET > 0.1, na.rm = TRUE),
    .groups    = "drop"
  ) %>%
  mutate(population = factor(population, levels = all_populations))

# =============================================================================
# 5. FST MATRIX from global Fst files
# =============================================================================

fst_files <- list.files(fst_dir,
                        pattern = "fst_pop_map_.*\\.fst\\.global\\.txt$",
                        full.names = TRUE)

fst_pairs <- lapply(fst_files, function(f) {
  pops <- sub(".*fst_pop_map_([A-Z]+)_pop_map_([A-Z]+)\\.fst\\.global\\.txt",
              "\\1_\\2", basename(f))
  pop1 <- strsplit(pops, "_")[[1]][1]
  pop2 <- strsplit(pops, "_")[[1]][2]
  lines <- readLines(f)
  last  <- trimws(tail(lines, 1))
  vals  <- as.numeric(strsplit(last, "\\s+")[[1]])
  tibble(pop1 = pop1, pop2 = pop2,
         Fst_unweighted = vals[1], Fst_weighted = vals[2])
}) %>%
  bind_rows()

# Build symmetric matrix
fst_pops   <- intersect(all_populations,
                        unique(c(fst_pairs$pop1, fst_pairs$pop2)))
fst_matrix <- matrix(NA, nrow = length(fst_pops),
                     ncol = length(fst_pops),
                     dimnames = list(fst_pops, fst_pops))
diag(fst_matrix) <- 0

for (i in seq_len(nrow(fst_pairs))) {
  p1 <- fst_pairs$pop1[i]; p2 <- fst_pairs$pop2[i]
  if (p1 %in% fst_pops && p2 %in% fst_pops) {
    fst_matrix[p1, p2] <- fst_pairs$Fst_weighted[i]
    fst_matrix[p2, p1] <- fst_pairs$Fst_weighted[i]
  }
}

cat("\nWeighted Fst matrix:\n")
print(round(fst_matrix, 4))
write.csv(round(fst_matrix, 5),
          file.path(out_dir, "fst_matrix_weighted.csv"))

# =============================================================================
# 6. COMBINE INTO FINAL TABLE
# =============================================================================

pop_stats <- ho_table %>%
  left_join(theta_list %>%
              select(population, theta_pi, theta_W, Tajima_D),
            by = "population") %>%
  left_join(sfs_stats %>%
              select(population, allelic_richness, private_alleles),
            by = "population") %>%
  left_join(fhet_summary %>%
              select(population, n_ind_fhet, mean_F_HET, n_inbred),
            by = "population") %>%
  arrange(match(as.character(population), all_populations))

cat("\nFull population statistics table:\n")
print(pop_stats)

write_tsv(pop_stats, file.path(out_dir, "population_statistics.tsv"))

# =============================================================================
# 7. FORMATTED DISPLAY TABLE
# =============================================================================

pop_stats_display <- pop_stats %>%
  transmute(
    Population              = as.character(population),
    N                       = n_ind_Ho,
    `Ho (per site)`         = round(Ho, 4),
    `θπ (per site)`         = formatC(theta_pi, format = "e", digits = 3),
    `θW (per site)`         = formatC(theta_W,  format = "e", digits = 3),
    `Tajima's D`            = round(Tajima_D, 3),
    `Allelic richness`      = round(allelic_richness, 2),
    `Private alleles`       = round(private_alleles, 1),
    `Mean F_HET`            = round(mean_F_HET, 4),
    `N inbred (F > 0.1)`   = n_inbred
  )

cat("\nFormatted population statistics table:\n")
print(pop_stats_display)

# HTML table with footnotes
pop_stats_display %>%
  kable("html",
        caption = paste0(
          "Population genetic statistics. ",
          "Ho = observed heterozygosity per site estimated from ",
          "genotype posterior probabilities in ANGSD across ~181k sites. ",
          "θπ and θW = per-site nucleotide diversity and Watterson's theta ",
          "estimated from the site frequency spectrum (ANGSD/realSFS). ",
          "Allelic richness rarefied to n=", n_ref / 2, " diploid individuals. ",
          "F_HET = inbreeding coefficient [1 - (H_obs/H_exp)]; ",
          "F > 0.1 indicates moderate inbreeding. ",
          "COI, DH, SP have very small sample sizes (n ≤ 3) and statistics ",
          "should be interpreted with caution."
        ),
        align = "lrrrrrrrrrr") %>%
  kable_styling(bootstrap_options = c("striped", "hover", "condensed"),
                full_width = FALSE) %>%
  row_spec(which(pop_stats_display$Population %in% c("COI","DH","SP")),
           color = "grey60", italic = TRUE) %>%
  footnote(
    general = paste0(
      "Ho: estimated from beagle genotype likelihoods (ANGSD). ",
      "θπ/θW/Tajima's D: estimated from SFS using realSFS/thetaStat. ",
      "Allelic richness: rarefied to smallest population (n=",
      n_ref / 2, " individuals). ",
      "Italic rows = small populations (n ≤ 3), interpret with caution."
    ),
    general_title = "Methods: "
  ) %>%
  save_kable(file.path(out_dir, "population_statistics_table.html"))

cat("HTML table saved.\n")

# =============================================================================
# 8. FST HEATMAP
# =============================================================================

fst_long <- as.data.frame(fst_matrix) %>%
  tibble::rownames_to_column("pop1") %>%
  pivot_longer(-pop1, names_to = "pop2", values_to = "Fst") %>%
  mutate(
    pop1 = factor(pop1, levels = all_populations),
    pop2 = factor(pop2, levels = rev(all_populations))
  )

p_fst <- ggplot(fst_long, aes(x = pop1, y = pop2, fill = Fst)) +
  geom_tile(colour = "white", linewidth = 0.5) +
  geom_text(
    aes(label = case_when(
      is.na(Fst)    ~ "",
      pop1 == pop2  ~ "—",
      TRUE          ~ sprintf("%.3f", Fst)
    )),
    size = 3.2
  ) +
  scale_fill_gradient(
    low      = "#E6F1FB",
    high     = "#0072B2",
    na.value = "grey95",
    name     = "Weighted\nFst"
  ) +
  labs(
    title    = "Pairwise weighted Fst between populations",
    subtitle = "Estimated from 2D-SFS using ANGSD/realSFS",
    x = NULL, y = NULL
  ) +
  theme_minimal(base_size = 12) +
  theme(axis.text  = element_text(face = "bold"),
        panel.grid = element_blank())

save_plot(p_fst, "fst_heatmap", width = 7, height = 6)

# =============================================================================
# 9. DIVERSITY COMPARISON PLOTS
# =============================================================================

# Filter to plot populations
plot_stats <- pop_stats %>%
  filter(population %in% plot_populations) %>%
  mutate(population = factor(population, levels = plot_populations))

# Ho bar plot
p_ho <- ggplot(plot_stats,
               aes(x = population, y = Ho, fill = population)) +
  geom_col(width = 0.6, alpha = 0.85) +
  scale_fill_manual(values = pop_colours) +
  labs(title    = "Observed heterozygosity (Ho) per population",
       subtitle = "Per-site Ho from genotype posterior probabilities (~181k sites)",
       x = NULL, y = "Ho (per site)") +
  theme_bw(base_size = 12) +
  theme(legend.position = "none")

# theta_pi bar plot
p_pi <- ggplot(plot_stats,
               aes(x = population, y = theta_pi, fill = population)) +
  geom_col(width = 0.6, alpha = 0.85) +
  scale_fill_manual(values = pop_colours) +
  labs(title    = "Nucleotide diversity (θπ) per population",
       subtitle = "Per-site θπ from site frequency spectrum (SFS)",
       x = NULL, y = "θπ (per site)") +
  theme_bw(base_size = 12) +
  theme(legend.position = "none")

# Tajima's D
p_tajima <- ggplot(plot_stats,
                   aes(x = population, y = Tajima_D, fill = population)) +
  geom_col(width = 0.6, alpha = 0.85) +
  geom_hline(yintercept = 0, linetype = "dashed",
             colour = "grey40", linewidth = 0.5) +
  scale_fill_manual(values = pop_colours) +
  labs(title    = "Tajima's D per population",
       subtitle = "Mean across windows; negative = population expansion or selection",
       x = NULL, y = "Tajima's D") +
  theme_bw(base_size = 12) +
  theme(legend.position = "none")

# Combined diversity panel
p_diversity <- (p_ho | p_pi | p_tajima) +
  plot_annotation(
    title = "Population diversity summary",
    theme = theme(plot.title = element_text(size = 14, face = "bold"))
  )

save_plot(p_diversity, "population_diversity_panel", width = 15, height = 5)
save_plot(p_ho,      "Ho_per_population",      width = 6, height = 5)
save_plot(p_pi,      "theta_pi_per_population", width = 6, height = 5)
save_plot(p_tajima,  "tajima_D_per_population", width = 6, height = 5)

cat("\n=== Population statistics complete ===\n")
cat("Output files in:", out_dir, "\n")
cat("  population_statistics.tsv\n")
cat("  population_statistics_table.html\n")
cat("  fst_matrix_weighted.csv\n")
cat("  fst_heatmap.pdf/.svg\n")
cat("  population_diversity_panel.pdf/.svg\n")
cat("  Ho_per_population.pdf/.svg\n")
cat("  theta_pi_per_population.pdf/.svg\n")
cat("  tajima_D_per_population.pdf/.svg\n")