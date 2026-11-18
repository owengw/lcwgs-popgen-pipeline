# =============================================================================
# NGSadmix admixture visualisation
# Features:
#   - CLUMPP-style component alignment across K values (Hungarian algorithm
#     via the clue package) so colours are consistent across K panels
#   - All 198 individual labels rotated 90 degrees on x-axis
#   - SI outliers and duplicate pair highlighted in red/brown
#   - Log-likelihood K selection plot with SD ribbon and all replicates
#   - Structure-style ancestry proportion plots
#   - Combined panel output
# All outputs as PDF and SVG
# =============================================================================

library(dplyr)
library(tidyr)
library(readr)
library(ggplot2)
library(patchwork)
library(stringr)
library(clue)       # Hungarian algorithm — install.packages("clue")

# =============================================================================
# USER SETTINGS
# =============================================================================

base_dir     <- file.path(Sys.getenv("USERPROFILE"),
                          "OneDrive - MMU/Documents/PhD Analysis/Bioinformatics/Physalia")
ngsadmix_dir <- file.path(base_dir, "admixture")
bam_list     <- file.path(base_dir, "pop_map_ALL.list")
meta_file    <- file.path(base_dir, "metadata.tsv")
out_dir      <- file.path(base_dir, "results", "admixture")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# Population order — consistent with PCA and diversity scripts
pop_order <- c("CAI", "NS", "SB", "SI", "COI", "SP")

# K values to plot structure bars for
K_TO_PLOT <- c(2, 3, 4, 5)

# Individuals to highlight (SI PCA outliers + duplicate pair)
FLAG_IDS <- c("1-SI75", "3-SI77", "8-SI82", "9-SI83")

# Component colours — assigned after alignment so consistent across K
component_colours <- c(
  "#F0E442",
  "#56B4E9",# comp1 — SI (yellow)
  "#E69F00",  # comp2 — CAI ancestry (orange) — SB/SP/DH also fall here
    # comp3 — NS (blue)
  "#CC79A7",  # comp4
  "#009E73",  # comp5
  "#0072B2",  # comp6
  "#D55E00",  # comp7
  "#999999",  # comp8
  "#44AA99",  # comp9
  "#882255"   # comp10
)

save_plot <- function(plot, name, width = 14, height = 6) {
  ggsave(file.path(out_dir, paste0(name, ".pdf")), plot,
         width = width, height = height)
  ggsave(file.path(out_dir, paste0(name, ".svg")), plot,
         width = width, height = height, device = "svg")
  cat("Saved:", name, "\n")
}

# =============================================================================
# 1. BUILD SAMPLE METADATA
# =============================================================================

cat("Loading sample metadata...\n")

list_lines <- read_lines(bam_list)

extract_id <- function(path) sub("_merged.*", "", sub(".*/", "", path))

samples <- tibble(
  row_idx   = seq_along(list_lines),
  bam_path  = list_lines,
  sample_id = extract_id(list_lines)
)

metadata <- read_tsv(meta_file, col_types = cols(), show_col_types = FALSE)

samples <- samples %>%
  mutate(sample_prefix = str_extract(sample_id, "^[^-]+-[A-Z]+")) %>%
  left_join(
    metadata %>%
      mutate(sample_prefix = str_extract(sample_id, "^[^-]+-[A-Z]+")) %>%
      select(sample_prefix, site, year),
    by = "sample_prefix"
  ) %>%
  mutate(
    site    = factor(site, levels = pop_order),
    is_flag = sample_id %in% FLAG_IDS,
    is_dup  = grepl("o$", str_extract(sample_id, "[A-Za-z0-9]+$")),
    label   = str_extract(sample_id, "[A-Z]+[0-9]+[a-z]*$")
  )

# Order individuals: by population then by sample_id within population
samples_ordered <- samples %>%
  arrange(site, sample_id) %>%
  mutate(plot_pos = row_number())

# Population midpoints and boundaries for panel annotations
pop_meta <- samples_ordered %>%
  group_by(site) %>%
  summarise(
    xmin = min(plot_pos) - 0.5,
    xmax = max(plot_pos) + 0.5,
    xmid = mean(plot_pos),
    .groups = "drop"
  )

# X-axis label colours: red = flagged outliers, brown = originals, grey = normal
label_colours <- samples_ordered %>%
  mutate(col = case_when(
    is_flag ~ "#A32D2D",
    is_dup  ~ "#854F0B",
    TRUE    ~ "grey40"
  )) %>%
  pull(col)

cat("Total individuals:", nrow(samples), "\n")
cat("Flagged (outliers):", sum(samples$is_flag), "\n")
cat("Original (non-amplified):", sum(samples$is_dup), "\n")

# =============================================================================
# 2. LOG-LIKELIHOOD K SELECTION
# =============================================================================

cat("\nLoading log-likelihoods...\n")

loglike_file <- file.path(ngsadmix_dir, "loglikelihoods.txt")
if (!file.exists(loglike_file)) stop("loglikelihoods.txt not found: ", loglike_file)

loglike <- read_lines(loglike_file) %>%
  tibble(raw = .) %>%
  mutate(
    K       = as.integer(str_extract(raw, "(?<=K=)\\d+")),
    rep     = as.integer(str_extract(raw, "(?<=Rep=)\\d+")),
    loglike = as.numeric(str_extract(raw, "(?<=loglike=)-?[0-9.]+"))
  ) %>%
  filter(!is.na(K), !is.na(loglike))

loglike_summary <- loglike %>%
  group_by(K) %>%
  summarise(
    n_reps  = n(),
    mean_ll = mean(loglike),
    sd_ll   = sd(loglike),
    max_ll  = max(loglike),
    min_ll  = min(loglike),
    .groups = "drop"
  )

best_reps <- loglike %>%
  group_by(K) %>%
  slice_max(loglike, n = 1, with_ties = FALSE) %>%
  ungroup()

# Last stable K = last K where SD < 10 (replicates converge consistently)
sd_stable <- loglike_summary %>%
  filter(sd_ll < 10) %>%
  slice_max(K, n = 1)

cat("Last stable K (SD < 10):", sd_stable$K, "\n")

p_ksel <- ggplot(loglike_summary, aes(x = K)) +
  geom_ribbon(aes(ymin = mean_ll - sd_ll,
                  ymax = mean_ll + sd_ll),
              fill = "#56B4E9", alpha = 0.2) +
  geom_jitter(data    = loglike,
              aes(y   = loglike),
              width   = 0.06, size = 1.2,
              alpha   = 0.4, colour = "grey60", shape = 16) +
  geom_line(aes(y = mean_ll), colour = "#185FA5", linewidth = 0.9) +
  geom_point(aes(y = mean_ll), colour = "#185FA5", size = 3.5) +
  geom_point(data = best_reps,
             aes(y = loglike), shape = 17,
             colour = "#A32D2D", size = 3) +
  geom_vline(xintercept = sd_stable$K, linetype = "dashed",
             colour = "#009E73", linewidth = 0.7) +
  annotate("text",
           x = sd_stable$K + 0.12,
           y = min(loglike_summary$mean_ll - loglike_summary$sd_ll) * 0.9999,
           label = paste0("Last stable K = ", sd_stable$K),
           hjust = 0, size = 3.2, colour = "#009E73") +
  scale_x_continuous(breaks = sort(unique(loglike$K))) +
  labs(
    title    = "NGSadmix K selection",
    subtitle = paste0(
      "Line = mean log-likelihood ± SD  |  ",
      "grey = all replicates  |  triangles = best replicate  |  ",
      "green dashed = last stable K (SD < 10)"
    ),
    x = "K (number of ancestral populations)",
    y = "Log-likelihood"
  ) +
  theme_bw(base_size = 12) +
  theme(panel.grid.minor = element_blank())

save_plot(p_ksel, "ngsadmix_K_selection", width = 10, height = 5)

# =============================================================================
# 3. LOAD Q MATRICES
# =============================================================================

load_qopt <- function(K_val, rep_num) {
  f <- file.path(
    ngsadmix_dir,
    sprintf("K%d_rep%d", K_val, rep_num),
    sprintf("ngsadmix_K%d_rep%d.qopt", K_val, rep_num)
  )
  if (!file.exists(f)) { cat("WARNING: Not found:", f, "\n"); return(NULL) }
  as.matrix(read_table(f, col_names = FALSE,
                       col_types = cols(.default = col_double())))
}

q_matrices <- list()
for (K_val in K_TO_PLOT) {
  best_rep <- filter(best_reps, K == K_val)$rep
  if (length(best_rep) == 0) next
  q_mat <- load_qopt(K_val, best_rep)
  if (!is.null(q_mat)) {
    q_matrices[[as.character(K_val)]] <- q_mat
    cat("Loaded K=", K_val, "rep=", best_rep,
        " dims:", nrow(q_mat), "x", ncol(q_mat), "\n")
  }
}

# =============================================================================
# 4. CLUMPP-STYLE ALIGNMENT (Hungarian algorithm via clue::solve_LSAP)
#
# Anchor: lowest K in K_TO_PLOT (most stable solution)
# For each higher K: find column permutation of q_mat that maximises
# correlation with the anchor matrix columns.
# Extra components (beyond anchor K) are appended after aligned ones.
# =============================================================================

cat("\nAligning components across K values...\n")

anchor_K   <- min(K_TO_PLOT)
anchor_mat <- q_matrices[[as.character(anchor_K)]]
K_anchor   <- ncol(anchor_mat)

aligned_matrices <- list()
aligned_matrices[[as.character(anchor_K)]] <- anchor_mat

for (K_val in K_TO_PLOT[K_TO_PLOT > anchor_K]) {
  q_mat  <- q_matrices[[as.character(K_val)]]
  if (is.null(q_mat)) next
  K_curr <- ncol(q_mat)
  
  # Cost matrix (K_curr x K_curr):
  # entry [i,j] = negative correlation between anchor col i and current col j
  # Rows beyond K_anchor are padded with 0 (no preference for unmatched)
  cost_mat <- matrix(0, nrow = K_curr, ncol = K_curr)
  for (i in seq_len(K_anchor)) {
    for (j in seq_len(K_curr)) {
      val <- cor(anchor_mat[, i], q_mat[, j])
      cost_mat[i, j] <- if (is.na(val)) 0 else -val
    }
  }
  
  cost_mat <- cost_mat - min(cost_mat)  # shift so minimum is 0
  cat("  Cost matrix NAs:", sum(is.na(cost_mat)), "\n")
  cat("  Cost matrix range:", min(cost_mat), "to", max(cost_mat), "\n")
  assignment      <- as.integer(solve_LSAP(cost_mat))
  aligned_mat     <- q_mat[, assignment]
  aligned_matrices[[as.character(K_val)]] <- aligned_mat
  
  cat("  K=", K_val, "permutation:", assignment, "\n")
}

# =============================================================================
# 5. STRUCTURE PLOT FUNCTION
# =============================================================================

plot_structure <- function(K_val, q_aligned) {
  K_curr     <- ncol(q_aligned)
  comp_names <- paste0("comp", seq_len(K_curr))
  colnames(q_aligned) <- comp_names
  
  q_df <- bind_cols(samples_ordered, as_tibble(q_aligned)) %>%
    pivot_longer(cols      = all_of(comp_names),
                 names_to  = "component",
                 values_to = "proportion") %>%
    mutate(component = factor(component, levels = comp_names))
  
  comp_cols <- setNames(component_colours[seq_len(K_curr)], comp_names)
  
  best_ll <- filter(best_reps, K == K_val)$loglike
  
  ggplot(q_df, aes(x = plot_pos, y = proportion, fill = component)) +
    geom_col(width = 1, position = "stack", linewidth = 0) +
    geom_vline(data        = pop_meta,
               aes(xintercept = xmax),
               colour      = "white", linewidth = 0.8,
               inherit.aes = FALSE) +
    geom_text(data        = pop_meta,
              aes(x = xmid, y = 1.08, label = site),
              inherit.aes = FALSE,
              size = 3.5, fontface = "bold", colour = "grey20") +
    scale_fill_manual(values = comp_cols, guide = "none") +
    scale_y_continuous(expand = c(0, 0),
                       limits = c(0, 1.14),
                       breaks = c(0, 0.25, 0.5, 0.75, 1.0)) +
    scale_x_continuous(
      expand = c(0, 0),
      breaks = samples_ordered$plot_pos,
      labels = samples_ordered$label
    ) +
    labs(
      title = sprintf("K = %d  (loglike = %.1f)", K_val, best_ll),
      x     = NULL,
      y     = "Ancestry\nproportion"
    ) +
    theme_bw(base_size = 10) +
    theme(
      axis.text.x  = element_text(angle = 90, hjust = 1, vjust = 0.5,
                                  size = 5.5, colour = label_colours),
      axis.ticks.x = element_line(linewidth = 0.2),
      panel.grid   = element_blank(),
      panel.border = element_rect(colour = "grey70"),
      plot.margin  = margin(5, 5, 5, 5),
      plot.title   = element_text(size = 10, face = "bold")
    )
}

# =============================================================================
# 6. GENERATE AND SAVE STRUCTURE PLOTS
# =============================================================================

structure_plots <- list()

for (K_val in K_TO_PLOT) {
  q_al <- aligned_matrices[[as.character(K_val)]]
  if (is.null(q_al)) next
  cat("Plotting K=", K_val, "\n")
  p <- plot_structure(K_val, q_al)
  structure_plots[[as.character(K_val)]] <- p
  save_plot(p, sprintf("ngsadmix_structure_K%d", K_val),
            width = 16, height = 3.8)
}

# =============================================================================
# 7. COMBINED STRUCTURE PANEL (all K values stacked)
# =============================================================================

if (length(structure_plots) > 0) {
  p_panel <- wrap_plots(structure_plots, ncol = 1) +
    plot_annotation(
      title   = "NGSadmix ancestry proportions",
      subtitle = sprintf("K = %s | best replicate per K | components aligned via Hungarian algorithm",
                         paste(K_TO_PLOT, collapse = ", ")),
      caption = paste0(
        "Red labels = PCA outliers (SI75, SI77, SI82, SI83)  |  ",
        "Brown labels = original non-amplified samples."
      ),
      theme = theme(
        plot.title    = element_text(size = 13, face = "bold"),
        plot.subtitle = element_text(size = 10),
        plot.caption  = element_text(size = 8, colour = "grey50")
      )
    )
  
  save_plot(p_panel, "ngsadmix_structure_combined",
            width = 16, height = 3.8 * length(structure_plots) + 1)
}

# =============================================================================
# 8. COMBINED K SELECTION + BEST STABLE K STRUCTURE
# =============================================================================

if (as.character(sd_stable$K) %in% names(structure_plots)) {
  p_final <- (p_ksel / structure_plots[[as.character(sd_stable$K)]]) +
    plot_layout(heights = c(1, 0.85)) +
    plot_annotation(
      title = sprintf("NGSadmix — K selection and best supported K=%d",
                      sd_stable$K),
      theme = theme(plot.title = element_text(size = 13, face = "bold"))
    )
  save_plot(p_final, "ngsadmix_combined_panel", width = 16, height = 10)
}

# =============================================================================
# 9. SUMMARIES
# =============================================================================

write_tsv(loglike_summary,
          file.path(out_dir, "ngsadmix_loglike_summary.tsv"))
write_tsv(best_reps %>% select(K, rep, loglike),
          file.path(out_dir, "ngsadmix_best_reps.tsv"))

cat("\nDone. Plots saved to:", out_dir, "\n")
cat("Last stable K:", sd_stable$K, "\n")
cat("Tip: also inspect K=", sd_stable$K + 1, "for biological interpretability.\n")

#K=1, 2, 3 have essentially zero variance (SD < 0.15), meaning all replicates converge to the same solution. 
#At K=4 the SD explodes to 4,146 and stays high through K=10. 
#This is the classic signal that K=3 is the last stable solution
#beyond that the algorithm is finding multiple different local optima
#the data don't strongly support additional structure. 
#K=4 and above are biologically plausible but statistically unstable.