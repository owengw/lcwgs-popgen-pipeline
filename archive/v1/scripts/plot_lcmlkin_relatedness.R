# =============================================================================
# lcMLkin relatedness visualisation
# Network plots + relationship summary table
# Uses igraph base plotting (no ggraph dependency)
#
# Input:  lcmlkin_all_populations.tsv
# Output: results/relatedness_lcmlkin/ (PDF + PNG per plot)
#
# lcMLkin columns:
#   Ind1, Ind2    — full BAM paths
#   Z0ag (K0)     — P(0 alleles IBD)
#   Z1ag (K1)     — P(1 allele IBD)
#   Z2ag (K2)     — P(2 alleles IBD)
#   PI_HATag (r)  — relatedness coefficient
#   nbSNP         — SNPs used
#   population    — added by collation
#   relationship  — classified by collation
#
# Relationship thresholds (r):
#   Duplicate:    r >= 0.45
#   1st degree:   K0 <= 0.10, r >= 0.40
#   Full sibling: K0 0.15-0.35, r >= 0.35
#   2nd degree:   K0 0.40-0.60
#   3rd degree:   K0 0.65-0.85
#   Unrelated:    K0 ~1.0
# =============================================================================

library(dplyr)
library(tidyr)
library(ggplot2)
library(igraph)
library(stringr)
library(patchwork)

# =============================================================================
# USER SETTINGS
# =============================================================================

base_dir    <- file.path(Sys.getenv("USERPROFILE"),
                         "OneDrive - MMU/Documents/PhD Analysis/Bioinformatics/Physalia")
lcmlkin_dir <- file.path(base_dir, "lcmlkin")
out_dir     <- file.path(base_dir, "results", "relatedness_lcmlkin")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

pop_order   <- c("CAI", "NS", "SB", "SI")
pop_colours <- c(CAI = "#E69F00", NS = "#56B4E9", SB = "#009E73", SI = "#F0E442")

rel_colours <- c(
  "Duplicate"    = "#A32D2D",
  "1st degree"   = "#D55E00",
  "Full sibling" = "#CC79A7",
  "2nd degree"   = "#E69F00",
  "3rd degree"   = "#56B4E9",
  "Unrelated"    = "grey85"
)

# =============================================================================
# 1. LOAD AND CLEAN DATA
# =============================================================================

cat("Loading lcMLkin results...\n")

lcmlkin_file <- file.path(lcmlkin_dir, "lcmlkin_all_populations.tsv")
if (!file.exists(lcmlkin_file)) {
  lcmlkin_file <- file.path(base_dir, "angsd_results", "lcmlkin",
                            "lcmlkin_all_populations.tsv")
}
stopifnot(file.exists(lcmlkin_file))

dat_raw <- read.table(lcmlkin_file, header = TRUE, sep = "\t",
                      stringsAsFactors = FALSE)

cat("Rows loaded:", nrow(dat_raw), "\n")

clean_id <- function(x) {
  x <- basename(x)
  x <- sub("_merged_pe_bt2.*$", "", x)
  x <- sub("_merged_se_bt2.*$", "", x)
  x
}

dat <- dat_raw %>%
  filter(PI_HATag != -9, Z0ag != -9) %>%
  mutate(
    Ind1_clean = clean_id(Ind1),
    Ind2_clean = clean_id(Ind2),
    Ind1_label = sub("^[0-9]+-", "", Ind1_clean),
    Ind2_label = sub("^[0-9]+-", "", Ind2_clean),
    population = factor(population, levels = pop_order),
    relationship = case_when(
      PI_HATag >= 0.45                                ~ "Duplicate",
      Z0ag <= 0.10 & PI_HATag >= 0.40                ~ "1st degree",
      Z0ag >= 0.15 & Z0ag <= 0.35 & PI_HATag >= 0.35 ~ "Full sibling",
      Z0ag >= 0.40 & Z0ag <= 0.60                    ~ "2nd degree",
      Z0ag >= 0.65 & Z0ag <= 0.85                    ~ "3rd degree",
      TRUE                                            ~ "Unrelated"
    ),
    relationship = factor(relationship,
                          levels = c("Duplicate", "1st degree", "Full sibling",
                                     "2nd degree", "3rd degree", "Unrelated"))
  )

cat("Valid pairs after filtering -9:", nrow(dat), "\n")

# =============================================================================
# 2. SUMMARY TABLE
# =============================================================================

cat("\n=== Relationship class summary ===\n")

rel_summary <- dat %>%
  count(population, relationship) %>%
  pivot_wider(names_from = relationship, values_from = n, values_fill = 0)

print(rel_summary)

write.table(rel_summary,
            file.path(out_dir, "lcmlkin_relationship_summary.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)

close_pairs <- dat %>%
  filter(relationship != "Unrelated") %>%
  select(population, Ind1_label, Ind2_label,
         Z0ag, Z1ag, Z2ag, PI_HATag, nbSNP, relationship) %>%
  arrange(population, relationship, desc(PI_HATag))

write.table(close_pairs,
            file.path(out_dir, "lcmlkin_close_pairs.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)

cat("\nClose pairs:\n")
print(close_pairs)

# =============================================================================
# 3. SUMMARY BAR CHART (ggplot2)
# =============================================================================

close_summary <- dat %>%
  filter(relationship != "Unrelated") %>%
  count(population, relationship)

p_bar <- ggplot(close_summary,
                aes(x = population, y = n, fill = relationship)) +
  geom_col(position = "dodge", width = 0.7,
           colour = "grey30", linewidth = 0.3) +
  geom_text(aes(label = n),
            position = position_dodge(width = 0.7),
            vjust = -0.4, size = 3.2, colour = "grey20") +
  scale_fill_manual(values = rel_colours, name = "Relationship",
                    drop = TRUE) +
  labs(
    title    = "Close pairwise relationships per population (lcMLkin v2.1)",
    subtitle = "Unrelated pairs excluded",
    x        = "Population",
    y        = "Number of pairs"
  ) +
  theme_bw(base_size = 12) +
  theme(panel.grid.minor = element_blank(),
        legend.position  = "bottom")

ggsave(file.path(out_dir, "lcmlkin_relationship_counts.pdf"),
       p_bar, width = 8, height = 5)
ggsave(file.path(out_dir, "lcmlkin_relationship_counts.png"),
       p_bar, width = 8, height = 5, dpi = 300)
cat("Saved: lcmlkin_relationship_counts\n")

# =============================================================================
# 4. NETWORK PLOTS — igraph base plotting, one per population
# Nodes = individuals with >= 1 close pair
# Edges = close relationships (r >= 0.125, i.e. >= 3rd degree)
# Edge colour = relationship class
# Edge width = PI_HATag (r)
# Node colour = population colour
# =============================================================================

plot_network_igraph <- function(pop, min_r = 0.125,
                                pdf_file = NULL, png_file = NULL) {

  d <- dat %>%
    filter(population == pop,
           relationship != "Unrelated",
           PI_HATag >= min_r)

  if (nrow(d) == 0) {
    cat("  No close pairs for", pop, "— skipping\n")
    return(invisible(NULL))
  }

  # All individuals in this population (for node set)
  all_inds <- dat %>%
    filter(population == pop) %>%
    select(Ind1_label, Ind2_label) %>%
    unlist() %>%
    unique() %>%
    sort()

  # Build edge list
  edges <- d %>%
    select(from = Ind1_label, to = Ind2_label,
           r = PI_HATag, relationship) %>%
    mutate(
      edge_col   = rel_colours[as.character(relationship)],
      edge_width = scales::rescale(r, to = c(1, 5))
    )

  # Individuals involved in at least one close pair (shown as labelled nodes)
  close_inds <- unique(c(edges$from, edges$to))

  # Build igraph object — all individuals as nodes
  g <- graph_from_data_frame(
    edges[, c("from", "to")],
    directed = FALSE,
    vertices = data.frame(name = all_inds)
  )

  # Node attributes
  V(g)$color       <- pop_colours[[pop]]
  V(g)$frame.color <- "grey30"
  V(g)$size        <- ifelse(V(g)$name %in% close_inds, 10, 5)
  V(g)$label       <- ifelse(V(g)$name %in% close_inds, V(g)$name, NA)
  V(g)$label.cex   <- 0.65
  V(g)$label.color <- "grey10"
  V(g)$label.dist  <- 0.6

  # Edge attributes
  E(g)$color <- edges$edge_col
  E(g)$width <- edges$edge_width

  set.seed(42)
  lay <- layout_with_fr(g)

  # Legend items
  rels_present <- unique(as.character(d$relationship))
  leg_cols     <- rel_colours[rels_present]

  subtitle <- paste0(
    "n = ", length(all_inds), " individuals  |  ",
    nrow(edges), " pairs \u2265 r ", min_r
  )

  plot_fn <- function() {
    par(mar = c(2, 1, 3, 1), bg = "white")
    plot(g,
         layout          = lay,
         vertex.shape    = "circle",
         edge.curved     = 0.15,
         main            = pop,
         sub             = subtitle)
    legend("bottomleft",
           legend = rels_present,
           col    = leg_cols,
           lwd    = 3,
           bty    = "n",
           cex    = 0.8,
           title  = "Relationship")
    mtext(side = 3,
          text = "lcMLkin v2.1  |  edge colour = relationship  |  edge width = r",
          cex  = 0.7, col = "grey50", line = 0.1)
  }

  if (!is.null(pdf_file)) {
    pdf(pdf_file, width = 10, height = 9)
    plot_fn()
    dev.off()
    cat("Saved:", basename(pdf_file), "\n")
  }

  if (!is.null(png_file)) {
    png(png_file, width = 10, height = 9, units = "in", res = 300)
    plot_fn()
    dev.off()
    cat("Saved:", basename(png_file), "\n")
  }

  invisible(g)
}

for (pop in pop_order) {
  cat("Building network:", pop, "\n")
  plot_network_igraph(
    pop,
    pdf_file = file.path(out_dir, paste0("lcmlkin_network_", pop, ".pdf")),
    png_file = file.path(out_dir, paste0("lcmlkin_network_", pop, ".png"))
  )
}

# =============================================================================
# 5. COMBINED NETWORK PANEL (2x2 grid, PDF only)
# =============================================================================

cat("Building combined network panel...\n")

pdf(file.path(out_dir, "lcmlkin_network_combined.pdf"),
    width = 18, height = 16)
par(mfrow = c(2, 2), mar = c(2, 1, 3, 1), bg = "white")

for (pop in pop_order) {

  d <- dat %>%
    filter(population == pop,
           relationship != "Unrelated",
           PI_HATag >= 0.125)

  all_inds <- dat %>%
    filter(population == pop) %>%
    select(Ind1_label, Ind2_label) %>%
    unlist() %>% unique() %>% sort()

  if (nrow(d) == 0) {
    plot.new()
    title(main = pop, sub = "No close pairs")
    next
  }

  edges <- d %>%
    select(from = Ind1_label, to = Ind2_label,
           r = PI_HATag, relationship) %>%
    mutate(
      edge_col   = rel_colours[as.character(relationship)],
      edge_width = scales::rescale(r, to = c(1, 5))
    )

  close_inds <- unique(c(edges$from, edges$to))

  g <- graph_from_data_frame(
    edges[, c("from", "to")],
    directed = FALSE,
    vertices = data.frame(name = all_inds)
  )

  V(g)$color       <- pop_colours[[pop]]
  V(g)$frame.color <- "grey30"
  V(g)$size        <- ifelse(V(g)$name %in% close_inds, 9, 4)
  V(g)$label       <- ifelse(V(g)$name %in% close_inds, V(g)$name, NA)
  V(g)$label.cex   <- 0.55
  V(g)$label.color <- "grey10"
  V(g)$label.dist  <- 0.5
  E(g)$color       <- edges$edge_col
  E(g)$width       <- edges$edge_width

  set.seed(42)
  lay <- layout_with_fr(g)

  plot(g,
       layout       = lay,
       vertex.shape = "circle",
       edge.curved  = 0.15,
       main         = pop,
       sub          = paste0("n=", length(all_inds),
                             "  |  ", nrow(edges), " close pairs"))

  rels_present <- unique(as.character(d$relationship))
  legend("bottomleft",
         legend = rels_present,
         col    = rel_colours[rels_present],
         lwd    = 3, bty = "n", cex = 0.7,
         title  = "Relationship")
}

dev.off()
cat("Saved: lcmlkin_network_combined.pdf\n")

cat("\nAll outputs saved to:", out_dir, "\n")
