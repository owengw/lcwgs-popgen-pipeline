# =============================================================================
# GONE Results — Contemporary Ne over ~200 generations
# Plestiodon longirostris - Bermuda lcWGS
# Method: LD decay (Santiago et al. 2020, GONE v1.0)
# Generation time: 3 years (conservative estimate for Plestiodon longirostris)
#
# NOTE ON EXCLUDED POPULATIONS:
# CAI (n=22) and SB (n=17): GONE returned Ne estimates of ~12,000 and ~18,000
#   respectively, which are biologically implausible given census size estimates
#   of ~50 individuals from mark-recapture surveys. Small sample sizes are
#   known to inflate GONE Ne estimates due to insufficient LD information.
#   These populations are excluded from plots but included in supplementary.
# COI (n=3): Insufficient sample size for reliable LD-based Ne estimation.
# DH (n=1) and SP (n=3): Not run.
# =============================================================================

library(tidyverse)
library(patchwork)

POP_COLOURS <- c(
  CAI = "#E69F00", NS = "#56B4E9", SB = "#009E73",
  SI  = "#F0E442", COI = "#0072B2"
)

# Populations with reliable estimates (sufficient n, consistent with census size)
POPS_RELIABLE <- c("NS", "SI")

# All populations run (for supplementary)
POPS_ALL <- c("CAI", "NS", "SB", "SI", "COI")

GEN_TIME <- 3   # generation time in years

# =============================================================================
# FILE PATHS
# =============================================================================

BASE    <- "C:/Users/oweng/OneDrive - MMU/Documents/PhD Analysis/Bioinformatics/Physalia"
NE_DIR  <- file.path(BASE, "ne_estimate")
OUT_DIR <- file.path(NE_DIR, "plots")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

# =============================================================================
# SECTION 1: LOAD GONE RESULTS
# =============================================================================

cat("Loading GONE results...\n")

parse_gone <- function(pop, ne_dir) {
  f <- file.path(ne_dir, pop, paste0("Output_Ne_", pop))
  if (!file.exists(f)) {
    cat(sprintf("  MISSING: %s\n", f))
    return(NULL)
  }
  
  dat <- read_tsv(f, skip = 2, col_names = c("generation", "Ne"),
                  show_col_types = FALSE) %>%
    mutate(generation = as.numeric(generation),
           Ne         = as.numeric(Ne)) %>%
    filter(!is.na(generation), !is.na(Ne),
           is.finite(Ne), Ne > 0, Ne < 1e8) %>%
    mutate(
      population = pop,
      years_ago  = generation * GEN_TIME,
      reliable   = pop %in% POPS_RELIABLE
    )
  
  cat(sprintf("  %-5s (n=%s): %d generations, Ne range %s - %s %s\n",
              pop,
              c(CAI="22", NS="50", SB="17", SI="79", COI="3")[pop],
              nrow(dat),
              scales::comma(round(min(dat$Ne))),
              scales::comma(round(max(dat$Ne))),
              ifelse(pop %in% POPS_RELIABLE, "[RELIABLE]", "[EXCLUDED - see note]")))
  dat
}

gone_all <- map(POPS_ALL, ~parse_gone(.x, NE_DIR)) %>%
  compact() %>%
  bind_rows()

gone_reliable <- gone_all %>% filter(population %in% POPS_RELIABLE)

# Current Ne = most recent generation
current_all <- gone_all %>%
  group_by(population) %>%
  slice_min(generation, n = 1) %>%
  ungroup()

current_reliable <- current_all %>%
  filter(population %in% POPS_RELIABLE)

cat("\n  Current Ne summary:\n")
current_all %>%
  arrange(match(population, POPS_ALL)) %>%
  rowwise() %>%
  group_walk(~cat(sprintf("    %-5s  Ne = %s  %s\n",
                          .x$population,
                          scales::comma(round(.x$Ne)),
                          ifelse(.x$population %in% POPS_RELIABLE,
                                 "", "(excluded — see note)"))))

# =============================================================================
# SECTION 2: MAIN PLOTS (NS and SI only)
# =============================================================================

cat("\nGenerating plots...\n")

# --- Plot 1: NS and SI — generations on x-axis ---
p_main <- gone_reliable %>%
  ggplot(aes(x = generation, y = Ne, colour = population)) +
  geom_line(linewidth = 1.2) +
  scale_x_reverse(labels = scales::label_comma(),
                  breaks = c(1, 50, 100, 150, 200)) +
  scale_y_log10(labels = scales::label_comma()) +
  scale_colour_manual(values = POP_COLOURS, name = "Population") +
  labs(
    x        = "Generations ago",
    y        = expression(N[e] ~ "(log scale)"),
    title    = "Recent demographic history — GONE",
    subtitle = sprintf(
      "NS (n=50) and SI (n=79); generation time = %d years\nCAI, SB, COI excluded due to insufficient sample size",
      GEN_TIME)
  ) +
  theme_classic(base_size = 13) +
  theme(plot.title    = element_text(face = "bold"),
        plot.subtitle = element_text(colour = "grey40", size = 10),
        legend.position = "right")

# --- Plot 2: NS and SI — years before present ---
p_years <- gone_reliable %>%
  ggplot(aes(x = years_ago, y = Ne, colour = population)) +
  geom_line(linewidth = 1.2) +
  scale_x_reverse(labels = scales::label_comma()) +
  scale_y_log10(labels = scales::label_comma()) +
  scale_colour_manual(values = POP_COLOURS, name = "Population") +
  labs(
    x        = "Years before present",
    y        = expression(N[e] ~ "(log scale)"),
    title    = "Recent demographic history — GONE",
    subtitle = sprintf("Generation time = %d years assumed", GEN_TIME)
  ) +
  theme_classic(base_size = 13) +
  theme(plot.title    = element_text(face = "bold"),
        plot.subtitle = element_text(colour = "grey40", size = 10),
        legend.position = "right")

# --- Plot 3: Current Ne bar chart (reliable only) ---
p_current <- current_reliable %>%
  mutate(population = factor(population, levels = POPS_RELIABLE)) %>%
  ggplot(aes(x = population, y = Ne, fill = population)) +
  geom_col(alpha = 0.85, width = 0.5) +
  geom_text(aes(label = scales::comma(round(Ne))),
            vjust = -0.4, size = 4.5) +
  scale_fill_manual(values = POP_COLOURS) +
  scale_y_continuous(labels = scales::label_comma(),
                     expand = expansion(mult = c(0, 0.2))) +
  labs(
    x        = "Population",
    y        = expression(N[e]),
    title    = "Current effective population size",
    subtitle = "GONE estimate — most recent generation"
  ) +
  theme_classic(base_size = 13) +
  theme(plot.title      = element_text(face = "bold"),
        plot.subtitle   = element_text(colour = "grey40", size = 10),
        legend.position = "none")

# --- Plot 4: Faceted NS and SI ---
p_facet <- gone_reliable %>%
  ggplot(aes(x = generation, y = Ne, colour = population)) +
  geom_line(linewidth = 1.2) +
  scale_x_reverse(labels = scales::label_comma()) +
  scale_y_log10(labels = scales::label_comma()) +
  scale_colour_manual(values = POP_COLOURS) +
  facet_wrap(~population, ncol = 2, scales = "free_y") +
  labs(
    x     = "Generations ago",
    y     = expression(N[e] ~ "(log scale)"),
    title = "Demographic history — NS and SI"
  ) +
  theme_classic(base_size = 12) +
  theme(
    plot.title       = element_text(face = "bold"),
    legend.position  = "none",
    strip.background = element_rect(fill = "grey90", colour = NA),
    strip.text       = element_text(face = "bold", size = 12)
  )

# --- Combined main panel ---
p_combined <- (p_main | p_current) / (p_years | p_facet) +
  plot_annotation(
    title    = "Contemporary Effective Population Size — Plestiodon longirostris",
    subtitle = "GONE v1.0 (Santiago et al. 2020); NS and SI only (see methods for exclusions)",
    theme    = theme(
      plot.title    = element_text(face = "bold", size = 15),
      plot.subtitle = element_text(colour = "grey40", size = 11)
    )
  )

# =============================================================================
# SECTION 3: SUPPLEMENTARY PLOT (all populations)
# =============================================================================

p_supp <- gone_all %>%
  mutate(
    line_type = ifelse(population %in% POPS_RELIABLE, "solid", "dashed"),
    population = factor(population, levels = POPS_ALL)
  ) %>%
  ggplot(aes(x = generation, y = Ne, colour = population,
             linetype = reliable)) +
  geom_line(linewidth = 1) +
  scale_x_reverse(labels = scales::label_comma(),
                  breaks = c(1, 50, 100, 150, 200)) +
  scale_y_log10(labels = scales::label_comma()) +
  scale_colour_manual(values = POP_COLOURS, name = "Population") +
  scale_linetype_manual(values = c("TRUE" = "solid", "FALSE" = "dashed"),
                        labels = c("TRUE" = "Reliable", "FALSE" = "Excluded"),
                        name   = "Estimate") +
  labs(
    x        = "Generations ago",
    y        = expression(N[e] ~ "(log scale)"),
    title    = "GONE Ne estimates — all populations (supplementary)",
    subtitle = "Dashed = excluded due to small sample size or census size inconsistency"
  ) +
  theme_classic(base_size = 13) +
  theme(plot.title    = element_text(face = "bold"),
        plot.subtitle = element_text(colour = "grey40", size = 10))

# =============================================================================
# SECTION 4: SAVE OUTPUTS
# =============================================================================

cat("Saving outputs...\n")

write_tsv(gone_all,      file.path(OUT_DIR, "gone_all_populations.tsv"))
write_tsv(gone_reliable, file.path(OUT_DIR, "gone_reliable_populations.tsv"))
write_tsv(current_all,   file.path(OUT_DIR, "gone_current_ne_all.tsv"))

plots <- list(
  combined      = list(p = p_combined, w = 16, h = 14),
  main          = list(p = p_main,     w = 10, h = 7),
  years         = list(p = p_years,    w = 10, h = 7),
  current_ne    = list(p = p_current,  w = 6,  h = 6),
  faceted       = list(p = p_facet,    w = 10, h = 6),
  supplementary = list(p = p_supp,     w = 12, h = 7)
)

walk(names(plots), function(nm) {
  ggsave(file.path(OUT_DIR, paste0("gone_", nm, ".pdf")),
         plots[[nm]]$p, width = plots[[nm]]$w, height = plots[[nm]]$h,
         device = cairo_pdf)
  ggsave(file.path(OUT_DIR, paste0("gone_", nm, ".svg")),
         plots[[nm]]$p, width = plots[[nm]]$w, height = plots[[nm]]$h)
  cat(sprintf("  Saved: gone_%s.pdf/svg\n", nm))
})

# =============================================================================
# SECTION 5: SUMMARY
# =============================================================================

cat("\n============================================\n")
cat("  GONE Ne Summary\n")
cat("============================================\n")
cat(sprintf("  Method:          GONE v1.0 (LD decay)\n"))
cat(sprintf("  Generation time: %d years\n", GEN_TIME))
cat(sprintf("  Generations:     %d\n", max(gone_all$generation)))
cat("\n  Reliable estimates (NS, SI):\n")
gone_reliable %>%
  group_by(population) %>%
  slice_min(generation, n = 1) %>%
  ungroup() %>%
  rowwise() %>%
  group_walk(~cat(sprintf("    %-5s  Ne = %s\n",
                          .x$population,
                          scales::comma(round(.x$Ne)))))
cat("\n  Excluded populations:\n")
cat("    CAI  Ne = ~12,000  (n=22; inconsistent with census size ~50)\n")
cat("    SB   Ne = ~18,000  (n=17; inconsistent with census size ~50)\n")
cat("    COI  Ne = ~94      (n=3;  insufficient sample size)\n")
cat(sprintf("\n  Outputs: %s\n", OUT_DIR))
cat("============================================\n")