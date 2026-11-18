# =============================================================================
# PCA from ANGSD covariance matrix (pop_map_ALL.covMat)
# LD-pruned SNPs, all 198 individuals — duplicates removed in R
# NaN values in covMat replaced with 0 (standard for ANGSD output)
# Produces: PC1v2, PC1v3, PC2v3, scree plot, duplicate-labelled plot
# All outputs as PDF and SVG
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

base_dir <- file.path(Sys.getenv("USERPROFILE"), 
                      "OneDrive - MMU/Documents/PhD Analysis/Bioinformatics/Physalia")
pca_dir   <- file.path(base_dir, "PCA")
meta_file <- file.path(base_dir, "metadata.tsv")
out_dir   <- file.path(base_dir, "results")
dir.create(out_dir, showWarnings = FALSE)

# Population display order and colours (Okabe-Ito colourblind-safe)
pop_order   <- c("CAI", "NS", "SB", "SI", "COI", "DH", "SP")
pop_colours <- c(
  CAI = "#E69F00", NS  = "#56B4E9", SB  = "#009E73", SI  = "#F0E442",
  COI = "#0072B2", DH  = "#D55E00", SP  = "#CC79A7"
)
pop_shapes <- c(
  CAI = 16, NS = 17, SB = 15, SI = 18,
  COI = 8,  DH = 10, SP = 12
)

save_plot <- function(plot, name, width = 8, height = 7) {
  ggsave(file.path(out_dir, paste0(name, ".pdf")), plot,
         width = width, height = height)
  ggsave(file.path(out_dir, paste0(name, ".svg")), plot,
         width = width, height = height, device = "svg")
  cat("Saved:", name, "\n")
}

# =============================================================================
# 1. BUILD SAMPLE METADATA
# =============================================================================

list_lines <- read_lines(file.path(pca_dir, "pop_map_ALL.list"))

extract_sample_id <- function(path) {
  sub("_.*", "", sub(".*/", "", path))
}

all_samples <- tibble(
  list_index = seq_along(list_lines),
  bam_path   = list_lines,
  sample_id  = extract_sample_id(list_lines),
  is_dup     = grepl("(aT|[^a]T|o)$", extract_sample_id(list_lines))
)

metadata <- read_tsv(meta_file, show_col_types = FALSE)

all_samples <- all_samples %>%
  left_join(metadata %>% select(sample_id, population, site, year),
            by = "sample_id")

cat("Total samples in list:", nrow(all_samples), "\n")
cat("Duplicates identified:", sum(all_samples$is_dup), "\n")
cat("Retained after removing duplicates:", sum(!all_samples$is_dup), "\n")

retained_idx     <- which(!all_samples$is_dup)
retained_samples <- all_samples %>% filter(!is_dup)

# =============================================================================
# 2. LOAD COVARIANCE MATRIX
# =============================================================================

cov_mat_raw <- as.matrix(read.table(
  file.path(pca_dir, "pcangsd_ALL.cov"),
  header = FALSE
))

cat("\nCovariance matrix dimensions:", nrow(cov_mat_raw), "x",
    ncol(cov_mat_raw), "\n")

# PCAngsd covariance matrices do not contain NaN values
# (unlike ANGSD IBS covMat); this check is retained as a safeguard
# individuals share no covered sites — ANGSD outputs NaN rather than 0)
n_nan <- sum(is.nan(cov_mat_raw))
cat("NaN values in raw matrix:", n_nan, "\n")

if (n_nan > 0) {
  cat("Replacing", n_nan, "NaN values with 0",
      "(no shared coverage between those pairs)\n")
  cov_mat_raw[is.nan(cov_mat_raw)] <- 0
}

# Subset to non-duplicate individuals
cov_mat <- cov_mat_raw[retained_idx, retained_idx]
cat("After removing duplicates:", nrow(cov_mat), "x", ncol(cov_mat), "\n")
cat("NAs remaining:", sum(is.na(cov_mat)), "\n")

# =============================================================================
# 3. EIGENDECOMPOSITION
# =============================================================================

cat("Running eigendecomposition...\n")
eig <- eigen(cov_mat, symmetric = TRUE)

# Variance explained per PC
var_explained <- eig$values / sum(eig$values) * 100
cum_var       <- cumsum(var_explained)

n_pcs <- min(10, nrow(cov_mat))

var_df <- data.frame(
  PC      = paste0("PC", 1:n_pcs),
  Var_pct = round(var_explained[1:n_pcs], 2),
  Cum_pct = round(cum_var[1:n_pcs], 2)
)
cat("\nVariance explained by first", n_pcs, "PCs:\n")
print(var_df)

# Build PCA scores data frame
pca_df <- as.data.frame(eig$vectors[, 1:n_pcs])
colnames(pca_df) <- paste0("PC", 1:n_pcs)

pca_df <- bind_cols(
  retained_samples %>% select(sample_id, population, site, year),
  pca_df
) %>%
  mutate(
    population = factor(population, levels = pop_order),
    label      = sample_id
  )

write_tsv(pca_df, file.path(out_dir, "pca_scores.tsv"))
cat("PCA scores saved.\n")

# =============================================================================
# 4. PLOT FUNCTION
# =============================================================================

plot_pca <- function(df, pcx, pcy, var_exp,
                     title_suffix    = "",
                     label_outliers  = TRUE,
                     outlier_sd      = 3,
                     add_ellipses    = FALSE) {
  
  xvar  <- round(var_exp[pcx], 1)
  yvar  <- round(var_exp[pcy], 1)
  x_col <- paste0("PC", pcx)
  y_col <- paste0("PC", pcy)
  
  df <- df %>%
    mutate(
      x_val   = .data[[x_col]],
      y_val   = .data[[y_col]],
      outlier = abs(scale(x_val)) > outlier_sd |
        abs(scale(y_val)) > outlier_sd
    )
  
  p <- ggplot(df, aes(x = x_val, y = y_val,
                      colour = population,
                      shape  = population)) +
    geom_hline(yintercept = 0, colour = "grey80",
               linewidth = 0.4, linetype = "dashed") +
    geom_vline(xintercept = 0, colour = "grey80",
               linewidth = 0.4, linetype = "dashed") +
    geom_point(size = 2.5, alpha = 0.85) +
    scale_colour_manual(values = pop_colours, name = "Population",
                        drop = FALSE) +
    scale_shape_manual(values  = pop_shapes,  name = "Population",
                       drop = FALSE) +
    labs(
      title    = paste0("PCA — PC", pcx, " vs PC", pcy, title_suffix),
      subtitle = paste0("PCAngsd covariance matrix, LD-pruned SNPs, n=",
                        nrow(df), " individuals"),
      x        = paste0("PC", pcx, " (", xvar, "%)"),
      y        = paste0("PC", pcy, " (", yvar, "%)")
    ) +
    theme_bw(base_size = 12) +
    theme(
      legend.position  = "right",
      panel.grid.minor = element_blank()
    )
  
  if (add_ellipses) {
    p <- p +
      stat_ellipse(level = 0.95, linewidth = 0.5,
                   linetype = "dashed", show.legend = FALSE)
  }
  
  if (label_outliers && any(df$outlier, na.rm = TRUE)) {
    p <- p +
      geom_text_repel(
        data          = filter(df, outlier),
        aes(label     = sub("^[0-9]+-", "", label)),
        size          = 2.5,
        colour        = "grey30",
        box.padding   = 0.4,
        max.overlaps  = 20
      )
  }
  
  p
}

# =============================================================================
# 5. SCREE PLOT
# =============================================================================

scree_df <- tibble(
  PC_num     = 1:n_pcs,
  Var_pct    = var_explained[1:n_pcs],
  Cumulative = cum_var[1:n_pcs]
)

p_scree <- ggplot(scree_df, aes(x = PC_num)) +
  geom_col(aes(y = Var_pct), fill = "steelblue",
           alpha = 0.8, width = 0.6) +
  geom_line(aes(y = Cumulative), colour = "coral",
            linewidth = 0.8) +
  geom_point(aes(y = Cumulative), colour = "coral", size = 2.5) +
  scale_x_continuous(breaks = 1:n_pcs,
                     labels = paste0("PC", 1:n_pcs)) +
  scale_y_continuous(
    name   = "Variance explained (%)",
    limits = c(0, max(scree_df$Cumulative) * 1.05)
  ) +
  labs(
    title    = "PCA scree plot",
    subtitle = "Bars = variance per PC; line = cumulative variance (%)",
    x        = NULL
  ) +
  theme_bw(base_size = 12) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

save_plot(p_scree, "pca_scree", width = 8, height = 5)

# =============================================================================
# 6. PC BIPLOTS — plain and with ellipses
# =============================================================================

for (ellipses in c(FALSE, TRUE)) {
  suffix <- if (ellipses) "_ellipses" else ""
  
  p_12 <- plot_pca(pca_df, 1, 2, var_explained,
                   add_ellipses = ellipses)
  p_13 <- plot_pca(pca_df, 1, 3, var_explained,
                   add_ellipses = ellipses)
  p_23 <- plot_pca(pca_df, 2, 3, var_explained,
                   add_ellipses = ellipses)
  
  save_plot(p_12, paste0("pca_PC1_PC2", suffix))
  save_plot(p_13, paste0("pca_PC1_PC3", suffix))
  save_plot(p_23, paste0("pca_PC2_PC3", suffix))
  
  p_combined <- (p_12 | p_13) / (p_23 | p_scree) +
    plot_layout(guides = "collect") +
    plot_annotation(
      title = paste0("PCA — ANGSD covariance matrix (LD-pruned SNPs)",
                     if (ellipses) " with 95% ellipses" else ""),
      theme = theme(
        plot.title      = element_text(size = 14, face = "bold"),
        legend.position = "right"
      )
    )
  
  save_plot(p_combined, paste0("pca_combined_panel", suffix),
            width = 16, height = 14)
}

# =============================================================================
# 7. PCA WITH DUPLICATES INCLUDED AND LABELLED
# =============================================================================

eig_full      <- eigen(cov_mat_raw, symmetric = TRUE)
var_exp_full  <- eig_full$values / sum(eig_full$values) * 100

pca_full <- as.data.frame(eig_full$vectors[, 1:n_pcs])
colnames(pca_full) <- paste0("PC", 1:n_pcs)

pca_full <- bind_cols(
  all_samples %>% select(sample_id, population, site, year, is_dup),
  pca_full
) %>%
  mutate(
    population = factor(population, levels = pop_order),
    label      = if_else(is_dup,
                         sub("^[0-9o]+-", "", sample_id),
                         NA_character_)
  )

xvar_full <- round(var_exp_full[1], 1)
yvar_full <- round(var_exp_full[2], 1)

p_dup <- ggplot(pca_full,
                aes(x = PC1, y = PC2,
                    colour = population,
                    shape  = population)) +
  geom_hline(yintercept = 0, colour = "grey80",
             linewidth = 0.4, linetype = "dashed") +
  geom_vline(xintercept = 0, colour = "grey80",
             linewidth = 0.4, linetype = "dashed") +
  # Non-duplicate points first
  geom_point(data   = filter(pca_full, !is_dup),
             size   = 2, alpha = 0.7) +
  # Duplicate points on top with black outline
  geom_point(data   = filter(pca_full, is_dup),
             size   = 4, alpha = 1, stroke = 1,
             shape  = 21,
             aes(fill = population),
             colour = "black") +
  geom_text_repel(
    data           = filter(pca_full, is_dup),
    aes(label      = label),
    size           = 2.8,
    colour         = "black",
    fontface       = "bold",
    box.padding    = 0.5,
    point.padding  = 0.3,
    max.overlaps   = Inf,
    segment.colour = "grey50",
    segment.size   = 0.3
  ) +
  scale_colour_manual(values = pop_colours, name = "Population",
                      drop = FALSE) +
  scale_fill_manual(values   = pop_colours, name = "Population",
                    drop = FALSE) +
  scale_shape_manual(values  = pop_shapes,  name = "Population",
                     drop = FALSE) +
  labs(
    title    = "PCA including duplicate samples (labelled)",
    subtitle = paste0(
      "Filled circles with black outline = duplicate/replicate samples  |  ",
      "n=198 individuals (including 15 duplicates)"
    ),
    x = paste0("PC1 (", xvar_full, "%)"),
    y = paste0("PC2 (", yvar_full, "%)")
  ) +
  theme_bw(base_size = 12) +
  theme(
    legend.position  = "right",
    panel.grid.minor = element_blank()
  )

save_plot(p_dup, "pca_with_duplicates_labelled", width = 10, height = 8)

# =============================================================================
# 8. SUMMARY
# =============================================================================

cat("\n=== PCA complete ===\n")
cat("Output files in:", out_dir, "\n")
cat("  pca_scores.tsv\n")
cat("  pca_scree.pdf/.svg\n")
cat("  pca_PC1_PC2.pdf/.svg\n")
cat("  pca_PC1_PC3.pdf/.svg\n")
cat("  pca_PC2_PC3.pdf/.svg\n")
cat("  pca_PC1_PC2_ellipses.pdf/.svg\n")
cat("  pca_PC1_PC3_ellipses.pdf/.svg\n")
cat("  pca_PC2_PC3_ellipses.pdf/.svg\n")
cat("  pca_combined_panel.pdf/.svg\n")
cat("  pca_combined_panel_ellipses.pdf/.svg\n")
cat("  pca_with_duplicates_labelled.pdf/.svg\n")
