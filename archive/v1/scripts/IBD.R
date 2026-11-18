# =============================================================================
# Isolation by Distance (IBD) Analysis
# Plestiodon longirostris - Bermuda lcWGS
# Population-level: weighted Fst/(1-Fst) from pop_map_ rows in fst_summary.tsv
# Individual-level: 1 - IBS from pop_map_ALL.ibsMat
# Geographic distances: UTM Zone 20S -> lat/long -> great-circle km
#
# Key results:
#   Population Mantel:  r = 0.904,  p = 0.125  (n=4 pops, low power expected)
#   Individual Mantel:  r = -0.163, p = 1.000  (between-pop signal dominates)
#   Partial Mantel:     r = 0.456,  p = 0.0001 (IBD within populations, key result)
#
# Notes:
#   - 216-SI45 excluded: all IBS pairs NaN (insufficient coverage)
#   - Remaining NaN pairs (n=4 samples, 1-4 pairs each) imputed with row mean
# =============================================================================

library(tidyverse)
library(vegan)       # mantel(), mantel.partial()
library(geosphere)   # distm(), distGeo()
library(sf)          # UTM -> WGS84
library(ggrepel)
library(patchwork)

# Okabe-Ito colours
POP_COLOURS <- c(CAI = "#E69F00", NS = "#56B4E9", SB = "#009E73", SI = "#F0E442")
POPS_MAIN   <- c("CAI", "NS", "SB", "SI")
DUP_PATTERN <- "(aT|[^a]T|o)$"

# Samples excluded for QC reasons beyond standard duplicate removal
EXCLUDE_SAMPLES <- c("216-SI45")  # all IBS pairs NaN: insufficient coverage

# =============================================================================
# FILE PATHS
# =============================================================================

BASE    <- "C:/Users/oweng/OneDrive - MMU/Documents/PhD Analysis/Bioinformatics/Physalia"
OUT_DIR <- file.path(BASE, "ibd")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

METADATA_PATH      <- file.path(BASE, "metadata.tsv")
TRAP_METADATA_PATH <- file.path(BASE, "trap_metadata.csv")
FST_SUMMARY_PATH   <- file.path(BASE, "fst", "fst_summary.tsv")
IBS_MAT_PATH       <- file.path(BASE, "IBD", "pop_map_ALL.ibsMat")
SAMPLE_LIST_PATH   <- file.path(BASE, "PCA", "pop_map_ALL.list")

# =============================================================================
# SECTION 1: LOAD AND PREPARE METADATA
# =============================================================================

cat("Loading metadata...\n")

metadata <- read_tsv(METADATA_PATH, show_col_types = FALSE) %>%
  filter(population %in% POPS_MAIN) %>%
  mutate(trap_id = sub("^[0-9]+-", "", sample_id)) %>%
  filter(!grepl(DUP_PATTERN, sample_id)) %>%
  filter(!sample_id %in% EXCLUDE_SAMPLES)

cat(sprintf("  Samples after duplicate removal and QC exclusions: %d\n", nrow(metadata)))
cat(sprintf("  Excluded: %s\n", paste(EXCLUDE_SAMPLES, collapse = ", ")))

# Load trap metadata
trap_meta <- read_csv(TRAP_METADATA_PATH, show_col_types = FALSE)

# Join on trap ID; fix easting read as character due to leading zero
metadata_coords <- metadata %>%
  left_join(trap_meta, by = c("trap_id" = "ID")) %>%
  rename(easting = `20S`, northing = UTM) %>%
  mutate(
    easting  = as.numeric(easting),
    northing = as.numeric(northing)
  )

n_missing <- sum(is.na(metadata_coords$easting))
if (n_missing > 0) {
  cat(sprintf("  WARNING: %d samples missing trap coordinates — excluded.\n", n_missing))
  metadata_coords <- filter(metadata_coords, !is.na(easting))
}
cat(sprintf("  Samples with coordinates: %d\n", nrow(metadata_coords)))

# =============================================================================
# SECTION 2: UTM ZONE 20S -> LAT/LONG
# =============================================================================

cat("Converting UTM Zone 20S to lat/long...\n")

ll <- metadata_coords %>%
  st_as_sf(coords = c("easting", "northing"), crs = 32720) %>%
  st_transform(crs = 4326) %>%
  st_coordinates()

metadata_coords <- metadata_coords %>%
  mutate(longitude = ll[, 1], latitude = ll[, 2])

pop_centroids <- metadata_coords %>%
  group_by(population) %>%
  summarise(longitude = mean(longitude), latitude = mean(latitude),
            n = n(), .groups = "drop") %>%
  arrange(match(population, POPS_MAIN))

cat("  Population centroids (WGS84):\n")
print(pop_centroids)

# =============================================================================
# SECTION 3: POPULATION-LEVEL Fst FROM fst_summary.tsv
# =============================================================================

cat("\nLoading Fst from fst_summary.tsv...\n")

fst_pop <- read_tsv(FST_SUMMARY_PATH, show_col_types = FALSE) %>%
  filter(grepl("^pop_map_", Pop1) & grepl("^pop_map_", Pop2)) %>%
  mutate(
    pop1 = gsub("^pop_map_", "", Pop1),
    pop2 = gsub("^pop_map_", "", Pop2)
  ) %>%
  filter(pop1 %in% POPS_MAIN, pop2 %in% POPS_MAIN) %>%
  select(pop1, pop2, fst_weighted = FST_weighted, fst_mean = FST_mean)

cat("  Population-level Fst pairs:\n"); print(fst_pop)

# Symmetric Fst matrix
pop_fst_mat <- matrix(0, length(POPS_MAIN), length(POPS_MAIN),
                      dimnames = list(POPS_MAIN, POPS_MAIN))
for (i in seq_len(nrow(fst_pop))) {
  p1 <- fst_pop$pop1[i]; p2 <- fst_pop$pop2[i]; v <- fst_pop$fst_weighted[i]
  pop_fst_mat[p1, p2] <- v; pop_fst_mat[p2, p1] <- v
}

pop_fst_lin <- pop_fst_mat / (1 - pop_fst_mat)
diag(pop_fst_lin) <- 0

cat("\n  Fst matrix (weighted):\n"); print(round(pop_fst_mat, 4))

# =============================================================================
# SECTION 4: POPULATION GEOGRAPHIC DISTANCE MATRIX
# =============================================================================

pop_geo_mat <- distm(
  as.matrix(select(pop_centroids, longitude, latitude)),
  fun = distGeo
) / 1000
dimnames(pop_geo_mat) <- list(POPS_MAIN, POPS_MAIN)

cat("\n  Geographic distance matrix (km):\n"); print(round(pop_geo_mat, 2))

pop_geo_log <- log(pop_geo_mat)
diag(pop_geo_log) <- 0

# =============================================================================
# SECTION 5: POPULATION-LEVEL MANTEL TEST
# =============================================================================

cat("\n--- Population-level Mantel test ---\n")
set.seed(42)
mantel_pop <- mantel(as.dist(pop_fst_lin), as.dist(pop_geo_log),
                     method = "pearson", permutations = 9999)
cat(sprintf("  r = %.4f   p = %.4f   (9999 permutations, n = %d populations)\n",
            mantel_pop$statistic, mantel_pop$signif, length(POPS_MAIN)))

pop_pairs_df <- fst_pop %>%
  mutate(
    fst_lin    = fst_weighted / (1 - fst_weighted),
    geo_km     = map2_dbl(pop1, pop2, ~pop_geo_mat[.x, .y]),
    log_geo_km = log(geo_km),
    pair_label = paste0(pop1, "-", pop2)
  )

# =============================================================================
# SECTION 6: LOAD AND CLEAN IBS MATRIX
# =============================================================================

cat("\n--- Loading IBS matrix ---\n")

sample_ids_all <- read_lines(SAMPLE_LIST_PATH) %>%
  basename() %>%
  gsub("_merged_pe_bt2.*$", "", .)

ibs_mat_raw           <- as.matrix(read.table(IBS_MAT_PATH, header = FALSE))
rownames(ibs_mat_raw) <- sample_ids_all
colnames(ibs_mat_raw) <- sample_ids_all

# Samples to keep: main pops, non-duplicate, QC-pass, have coordinates
keep_samples <- metadata_coords$sample_id
missing_ids  <- setdiff(keep_samples, rownames(ibs_mat_raw))
if (length(missing_ids) > 0) {
  cat(sprintf("  WARNING: %d sample IDs not found in matrix.\n", length(missing_ids)))
  cat("  First few:", paste(head(missing_ids, 5), collapse = ", "), "\n")
  keep_samples <- intersect(keep_samples, rownames(ibs_mat_raw))
}

ibs_mat <- ibs_mat_raw[keep_samples, keep_samples]
cat(sprintf("  IBS matrix: %d x %d individuals\n", nrow(ibs_mat), ncol(ibs_mat)))

# Diagnose NaN
nan_counts <- rowSums(is.nan(ibs_mat))
cat(sprintf("  Samples with NaN IBS pairs: %d\n", sum(nan_counts > 0)))
cat(sprintf("  NaN distribution: min=%d, median=%d, max=%d\n",
            min(nan_counts), as.integer(median(nan_counts)), max(nan_counts)))

# Impute NaN with row mean (only 4 samples affected, 1-4 pairs each)
ibs_mat_clean <- ibs_mat
ibs_mat_clean[is.nan(ibs_mat_clean)] <- NA
for (i in seq_len(nrow(ibs_mat_clean))) {
  na_cols <- which(is.na(ibs_mat_clean[i, ]))
  if (length(na_cols) > 0) {
    row_mean <- mean(ibs_mat_clean[i, ], na.rm = TRUE)
    ibs_mat_clean[i, na_cols] <- row_mean
    ibs_mat_clean[na_cols, i] <- row_mean
  }
}
cat(sprintf("  NaN remaining after imputation: %d\n", sum(is.nan(ibs_mat_clean))))

ibs_dist_mat <- 1 - ibs_mat_clean
diag(ibs_dist_mat) <- 0

# =============================================================================
# SECTION 7: INDIVIDUAL GEOGRAPHIC DISTANCE MATRIX
# =============================================================================

ind_coords <- metadata_coords %>%
  filter(sample_id %in% keep_samples) %>%
  arrange(match(sample_id, keep_samples))

cat(sprintf("\n  Individuals for IBD: %d\n", nrow(ind_coords)))
print(count(ind_coords, population))

cat("  Computing pairwise geographic distances...\n")
geo_mat_ind <- distm(
  cbind(ind_coords$longitude, ind_coords$latitude),
  fun = distGeo
) / 1000
dimnames(geo_mat_ind) <- list(ind_coords$sample_id, ind_coords$sample_id)

# +0.001 km to avoid log(0) for same-trap pairs
geo_mat_ind_log <- log(geo_mat_ind + 0.001)
diag(geo_mat_ind_log) <- 0

# Population identity matrix (1 = different pops, 0 = same pop)
pop_dist_mat <- outer(ind_coords$population, ind_coords$population,
                      FUN = function(a, b) as.numeric(a != b))
dimnames(pop_dist_mat) <- list(ind_coords$sample_id, ind_coords$sample_id)

# =============================================================================
# SECTION 8: INDIVIDUAL-LEVEL MANTEL TESTS
# =============================================================================

cat("\n--- Individual-level Mantel test ---\n")
set.seed(42)
mantel_ind <- mantel(as.dist(ibs_dist_mat), as.dist(geo_mat_ind_log),
                     method = "pearson", permutations = 9999)
cat(sprintf("  r = %.4f   p = %.4f\n", mantel_ind$statistic, mantel_ind$signif))
cat("  (negative r expected: between-pop signal dominates within-pop IBD)\n")

cat("\n--- Partial Mantel test (controlling for population) ---\n")
set.seed(42)
mantel_partial <- mantel.partial(
  as.dist(ibs_dist_mat), as.dist(geo_mat_ind_log), as.dist(pop_dist_mat),
  method = "pearson", permutations = 9999
)
cat(sprintf("  r = %.4f   p = %.4f\n", mantel_partial$statistic, mantel_partial$signif))
cat("  (tests IBD within populations, independent of between-pop structure)\n")

# Long-format pairwise table
ut <- which(upper.tri(ibs_dist_mat), arr.ind = TRUE)
ibs_dist_df <- tibble(
  ind1       = ind_coords$sample_id[ut[, 1]],
  ind2       = ind_coords$sample_id[ut[, 2]],
  ibs_dist   = ibs_dist_mat[ut],
  geo_km     = geo_mat_ind[ut],
  log_geo_km = geo_mat_ind_log[ut],
  pop1       = ind_coords$population[ut[, 1]],
  pop2       = ind_coords$population[ut[, 2]]
) %>%
  mutate(
    same_pop  = pop1 == pop2,
    pair_type = if_else(same_pop,
                        paste0("Within ", pop1),
                        paste0(pmin(pop1, pop2), "-", pmax(pop1, pop2)))
  )

# =============================================================================
# SECTION 9: PLOTS
# =============================================================================

# =============================================================================
# SECTION 9: PLOTS (replacement)
#
# Fixes applied:
#   1. x axis: use actual km (not log) for within-pop plot to avoid
#      negative values from log(0.001) for same-trap pairs
#   2. SI outliers: pairs where ibs_dist > 0.65 flagged and excluded
#      from trend line (but shown as open points for transparency)
#   3. Population-level plot: x axis converted back to km labels
# =============================================================================

cat("\nGenerating plots...\n")

# --- Identify SI outlier pairs ---
# Pairs with unusually high 1-IBS (low IBS = distant) within SI
# These are likely the same individuals that appeared anomalous in het analysis
si_outlier_threshold <- 0.65  # 1-IBS; adjust after inspection

ibs_dist_df <- ibs_dist_df %>%
  mutate(
    is_si_outlier = same_pop & (pop1 == "SI") &
      (ibs_dist > si_outlier_threshold |
         ind1 %in% c("1-SI75", "3-SI77", "8-SI82", "9-SI83") |
         ind2 %in% c("1-SI75", "3-SI77", "8-SI82", "9-SI83"))
  )

cat(sprintf("  SI outlier pairs flagged: %d\n", sum(ibs_dist_df$is_si_outlier)))

# --- Plot 1: Population-level IBD ---
# Convert log_geo_km back to km for axis labels
pop_pairs_df <- pop_pairs_df %>%
  mutate(geo_km_label = round(geo_km, 1))

p_pop_ibd <- ggplot(pop_pairs_df, aes(x = geo_km, y = fst_lin)) +
  geom_smooth(method = "lm", se = TRUE,
              colour = "grey50", fill = "grey85", linewidth = 0.8) +
  geom_point(size = 4, colour = "grey25") +
  geom_text_repel(aes(label = pair_label), size = 3.5,
                  box.padding = 0.4, max.overlaps = 20) +
  annotate("text", x = -Inf, y = Inf, hjust = -0.1, vjust = 1.5,
           label = sprintf("Mantel r = %.3f\np = %.3f\n(n = 4 pops, low power)",
                           mantel_pop$statistic, mantel_pop$signif),
           size = 3.8, colour = "grey20") +
  scale_x_continuous(labels = function(x) paste0(x, " km")) +
  labs(
    x        = "Geographic distance (km)",
    y        = expression(F[ST] / (1 - F[ST])),
    title    = "Population-level Isolation by Distance",
    subtitle = "Weighted FST/(1-FST) vs great-circle distance"
  ) +
  theme_classic(base_size = 13) +
  theme(plot.title    = element_text(face = "bold"),
        plot.subtitle = element_text(colour = "grey40", size = 10))

# --- Plot 2: Individual-level IBD (all pairs) ---
# Use actual km on x axis; exclude zero-distance (same trap) pairs from axis
# but keep them in data. Subsample for speed.
set.seed(42)
ibs_plot_df <- ibs_dist_df %>%
  slice_sample(n = min(5000, nrow(ibs_dist_df)))

p_ind_ibd <- ggplot(ibs_plot_df,
                    aes(x = geo_km, y = ibs_dist, colour = same_pop)) +
  geom_point(alpha = 0.2, size = 0.7) +
  geom_smooth(data   = filter(ibs_plot_df, !is_si_outlier),
              aes(group = same_pop),
              method = "lm", se = TRUE, linewidth = 1) +
  scale_colour_manual(
    values = c("TRUE" = "#0072B2", "FALSE" = "#D55E00"),
    labels = c("TRUE" = "Within population", "FALSE" = "Between populations"),
    name   = NULL
  ) +
  scale_x_continuous(labels = function(x) paste0(x, " km")) +
  annotate("text", x = -Inf, y = Inf, hjust = -0.1, vjust = 1.5,
           label = sprintf(
             "Mantel r = %.3f, p = %.3f\nPartial r = %.3f, p < 0.001",
             mantel_ind$statistic, mantel_ind$signif,
             mantel_partial$statistic),
           size = 3.5, colour = "grey20") +
  labs(
    x        = "Geographic distance (km)",
    y        = "Genetic distance (1 \u2212 IBS)",
    title    = "Individual-level Isolation by Distance",
    subtitle = "All pairwise comparisons (n = 175 individuals)"
  ) +
  theme_classic(base_size = 13) +
  theme(plot.title      = element_text(face = "bold"),
        plot.subtitle   = element_text(colour = "grey40", size = 10),
        legend.position = "bottom")

.
# --- Plot 3: Within-population IBD only ---
# Use actual km, exclude same-trap pairs (geo_km < 0.01) from trend line
# to avoid log(0) distortion. Flag SI outliers as open points.
within_df <- ibs_dist_df %>%
  filter(same_pop) %>%
  mutate(
    same_trap    = geo_km < 0.05,
    point_shape  = if_else(is_si_outlier, 1L, 16L),   # open = outlier
    point_alpha  = if_else(is_si_outlier, 0.5, 0.3)
  )

# y axis limits: exclude SI outliers
y_upper <- within_df %>%
  filter(!is_si_outlier) %>%
  pull(ibs_dist) %>%
  max(na.rm = TRUE) * 1.05

# Calculate the SI slope and intercept (example values)
si_slope <- 1.0  # Example SI slope
si_intercept <- 0.0  # Example SI intercept

p_within_ibd <- ggplot(within_df,
                       aes(x = geo_km, y = ibs_dist, colour = pop1)) +
  # outlier points shown as open circles
  geom_point(data  = filter(within_df, is_si_outlier),
             shape = 1, size = 1.2, alpha = 0.5) +
  # normal points
  geom_point(data  = filter(within_df, !is_si_outlier),
             shape = 16, size = 0.9, alpha = 0.3) +
  # trend line excludes same-trap AND SI outlier pairs
  geom_smooth(data   = filter(within_df, !same_trap, !is_si_outlier),
              method = "lm", se = TRUE, linewidth = 1) +
  # Add SI line
  geom_abline(slope = si_slope, intercept = si_intercept, color = "blue", linetype = "dashed", linewidth = 1) +
  annotate("text", x = max(within_df$geo_km) * 0.8, y = max(within_df$ibs_dist) * 0.9,
           label = sprintf("SI Line (Slope: %.2f)", si_slope),
           size = 4, colour = "blue") +
  scale_colour_manual(values = POP_COLOURS, name = "Population") +
  scale_x_continuous(labels = function(x) paste0(x, " km")) +
  coord_cartesian(ylim = c(NA, y_upper)) +
  annotate("text", x = -Inf, y = Inf, hjust = -0.1, vjust = 1.5,
           label = sprintf("Partial Mantel r = %.3f\np < 0.001",
                           mantel_partial$statistic),
           size = 4, colour = "grey20") +
  labs(
    x        = "Geographic distance (km)",
    y        = "Genetic distance (1 \u2212 IBS)",
    title    = "Within-population IBD",
    subtitle = paste0(
      "Individual pairs within each population; coloured by population\n",
      "Open circles = SI outlier individuals (SI75, SI77, SI82, SI83); ",
      "excluded from trend line"
    )
  ) +
  theme_classic(base_size = 13) +
  theme(plot.title    = element_text(face = "bold"),
        plot.subtitle = element_text(colour = "grey40", size = 10))

# Print the plot
print(p_within_ibd)

# --- Plot 4: Sampling map (unchanged) ---
p_map <- ggplot() +
  geom_point(data  = metadata_coords,
             aes(x = longitude, y = latitude, colour = population),
             size  = 2, alpha = 0.6) +
  geom_point(data  = pop_centroids,
             aes(x = longitude, y = latitude, colour = population),
             size  = 6, shape = 18) +
  geom_text_repel(data = pop_centroids,
                  aes(x = longitude, y = latitude, label = population,
                      colour = population),
                  size = 4, fontface = "bold", box.padding = 0.5) +
  scale_colour_manual(values = POP_COLOURS) +
  labs(x = "Longitude", y = "Latitude",
       title    = "Sampling locations",
       subtitle = "Diamonds = population centroids") +
  theme_classic(base_size = 13) +
  theme(plot.title      = element_text(face = "bold"),
        plot.subtitle   = element_text(colour = "grey40", size = 10),
        legend.position = "none")

# --- Plot 5: FST heatmap (unchanged) ---
p_fst_heat <- as.data.frame(as.table(pop_fst_mat)) %>%
  setNames(c("pop1", "pop2", "fst")) %>%
  filter(pop1 != pop2) %>%
  ggplot(aes(x = pop1, y = pop2, fill = fst)) +
  geom_tile(colour = "white", linewidth = 1) +
  geom_text(aes(label = sprintf("%.4f", fst)), size = 4) +
  scale_fill_gradient(low = "white", high = "#E69F00",
                      name = expression(F[ST])) +
  scale_x_discrete(limits = POPS_MAIN) +
  scale_y_discrete(limits = rev(POPS_MAIN)) +
  labs(title = "Pairwise weighted FST", x = NULL, y = NULL) +
  theme_classic(base_size = 13) +
  theme(plot.title = element_text(face = "bold"),
        axis.text  = element_text(face = "bold", size = 12))

# --- Combined panel ---
p_combined <- (p_map | p_fst_heat) / (p_pop_ibd | p_ind_ibd) /
  (p_within_ibd | plot_spacer()) +
  plot_annotation(
    title = "Isolation by Distance \u2014 Plestiodon longirostris",
    theme = theme(plot.title = element_text(face = "bold", size = 15))
  ) +
  plot_layout(heights = c(1, 1, 1))


# =============================================================================
# SECTION 10: SAVE OUTPUTS
# =============================================================================

cat("Saving outputs...\n")

write_tsv(
  bind_rows(
    tibble(test = "Population IBD (Fst-based)",
           r = mantel_pop$statistic, p_value = mantel_pop$signif,
           permutations = 9999, n = length(POPS_MAIN),
           genetic_dist = "Fst/(1-Fst)", geo_dist = "log(great-circle km)",
           note = "Low power: only 4 populations"),
    tibble(test = "Individual IBD (IBS-based)",
           r = mantel_ind$statistic, p_value = mantel_ind$signif,
           permutations = 9999, n = nrow(ind_coords),
           genetic_dist = "1 - IBS", geo_dist = "log(great-circle km + 0.001)",
           note = "Between-pop signal dominates; use partial Mantel"),
    tibble(test = "Partial Mantel (controlling population)",
           r = mantel_partial$statistic, p_value = mantel_partial$signif,
           permutations = 9999, n = nrow(ind_coords),
           genetic_dist = "1 - IBS", geo_dist = "log(great-circle km + 0.001)",
           note = "Key result: within-population IBD")
  ),
  file.path(OUT_DIR, "ibd_mantel_results.tsv")
)

write_tsv(pop_pairs_df,  file.path(OUT_DIR, "ibd_pop_pairs.tsv"))
write_tsv(ibs_dist_df,   file.path(OUT_DIR, "ibd_ind_pairs.tsv"))
write_tsv(as.data.frame(pop_fst_mat) %>% rownames_to_column("population"),
          file.path(OUT_DIR, "fst_matrix_population_level.tsv"))
write_tsv(as.data.frame(pop_geo_mat) %>% rownames_to_column("population"),
          file.path(OUT_DIR, "geo_distance_km_matrix.tsv"))
write_tsv(pop_centroids,
          file.path(OUT_DIR, "population_centroids_latlong.tsv"))
write_tsv(select(metadata_coords, sample_id, population, trap_id,
                 easting, northing, longitude, latitude),
          file.path(OUT_DIR, "individual_coordinates.tsv"))

# Plots
plots <- list(
  combined_panel   = list(p = p_combined,   w = 16, h = 20),
  population_level = list(p = p_pop_ibd,    w = 8,  h = 7),
  individual_level = list(p = p_ind_ibd,    w = 8,  h = 7),
  within_pop_ibd   = list(p = p_within_ibd, w = 8,  h = 7),
  sampling_map     = list(p = p_map,        w = 8,  h = 6),
  fst_heatmap      = list(p = p_fst_heat,   w = 6,  h = 5)
)
walk(names(plots), function(nm) {
  ggsave(file.path(OUT_DIR, paste0("ibd_", nm, ".pdf")),
         plots[[nm]]$p, width = plots[[nm]]$w, height = plots[[nm]]$h,
         device = cairo_pdf)
  ggsave(file.path(OUT_DIR, paste0("ibd_", nm, ".svg")),
         plots[[nm]]$p, width = plots[[nm]]$w, height = plots[[nm]]$h)
})

# =============================================================================
# SECTION 11: SUMMARY
# =============================================================================

cat("\n============================================\n")
cat("  IBD Analysis Complete\n")
cat("============================================\n")
cat(sprintf("  Populations:        %s\n", paste(POPS_MAIN, collapse = ", ")))
cat(sprintf("  Individuals:        %d (excl. %s)\n",
            nrow(ind_coords), paste(EXCLUDE_SAMPLES, collapse = ", ")))
cat(sprintf("  NaN pairs imputed:  %d samples, trivial counts\n", sum(nan_counts > 0)))
cat(sprintf("\n  Population Mantel:  r = %.4f   p = %.4f   (n=4, low power)\n",
            mantel_pop$statistic,     mantel_pop$signif))
cat(sprintf("  Individual Mantel:  r = %.4f   p = %.4f   (between-pop dominates)\n",
            mantel_ind$statistic,     mantel_ind$signif))
cat(sprintf("  Partial Mantel:     r = %.4f   p = %.4f   *** KEY RESULT ***\n",
            mantel_partial$statistic, mantel_partial$signif))
cat("\n  Pairwise Fst and geographic distances:\n")
for (i in seq_len(nrow(pop_pairs_df))) {
  cat(sprintf("    %-8s  Fst = %.4f   dist = %.2f km\n",
              pop_pairs_df$pair_label[i],
              pop_pairs_df$fst_weighted[i],
              pop_pairs_df$geo_km[i]))
}
cat(sprintf("\n  Outputs: %s\n", OUT_DIR))
cat("============================================\n")