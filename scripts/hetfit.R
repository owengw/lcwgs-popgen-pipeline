# =============================================================================
# Heterozygosity-Fitness Regression
# Plestiodon longirostris - Bermuda lcWGS
#
# Three modelling approaches compared:
#   1. OLS on pre-calculated SMI with robust standard errors (HC1)
#   2. OLS on log(SMI)
#   3. Beta regression on SMI scaled to (0,1)
#
# Fixes applied vs previous version:
#   - SMI used directly from body_stats.csv (pre-calculated, bsma=1.044)
#     rather than recalculated — avoids bsma discrepancy
#   - Bounded SMI derived by scaling pre-calculated SMI to (0,1)
#   - HC1 robust SEs used instead of HC3 — HC3 undefined when hat=1
#     (Daniel's Head, n=1 in its site group, has perfect leverage)
#   - Singleton sites (n=1: Daniel's Head) excluded from site models
#     to avoid perfect collinearity; retained in het-only models
#   - Beta regression interaction model dropped (overparameterised)
#   - NaN weight handled: individuals with missing weight excluded
#     with diagnostic output showing which IDs are affected
# =============================================================================

library(dplyr)
library(tidyr)
library(ggplot2)
library(patchwork)
library(sandwich)
library(lmtest)
library(betareg)
# emmeans and multcomp loaded later in post-hoc section (auto-installed if needed)

# Resolve namespace conflicts — betareg/multcomp can mask dplyr functions
library(conflicted)
conflicts_prefer(dplyr::select)
conflicts_prefer(dplyr::filter)
conflicts_prefer(dplyr::mutate)
conflicts_prefer(dplyr::summarise)

# =============================================================================
# PATHS
# =============================================================================

base_dir <- file.path(Sys.getenv("USERPROFILE"),
                      "OneDrive - MMU/Documents/PhD Analysis/Bioinformatics/Physalia")
het_dir  <- file.path(base_dir, "heterozygosity")
fit_dir  <- file.path(base_dir, "hetfit")
out_dir  <- file.path(base_dir, "hetfit", "results")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

pop_colours <- c(
  "Castle Island"      = "#E69F00",
  "Nonsuch Island"     = "#56B4E9",
  "Sinky Bay"          = "#009E73",
  "Southampton Island" = "#F0E442",
  "Coopers Island"     = "#0072B2",
  "Spittal Pond"       = "#CC79A7",
  "Daniel's Head"      = "#999999"
)

# =============================================================================
# 1. LOAD HETEROZYGOSITY
# =============================================================================

cat("Loading .het files...\n")
het_files <- list.files(het_dir, pattern = "\\.het$", full.names = TRUE)
cat("Found", length(het_files), ".het files\n")

het_list <- lapply(het_files, function(f) {
  tryCatch(read.table(f, header = TRUE, sep = "\t", stringsAsFactors = FALSE),
           error = function(e) NULL)
})
het_raw <- do.call(rbind, het_list[!sapply(het_list, is.null)])

het_raw$ID <- sapply(het_raw$sample_id, function(x) {
  x <- sub("_.*", "", x)
  parts <- strsplit(x, "-")[[1]]
  if (length(parts) >= 2) parts[2] else x
})

het_clean <- het_raw %>%
  select(ID, heterozygosity, n_sites_total) %>%
  filter(n_sites_total > 10000, heterozygosity > 0) %>%
  group_by(ID) %>%
  slice_max(n_sites_total, n = 1, with_ties = FALSE) %>%
  ungroup()

cat("Individuals after filtering:", nrow(het_clean), "\n")

# =============================================================================
# 2. LOAD BODY STATS
# Use pre-calculated SMI directly (bsma = 1.044, Peig & Green 2009)
# =============================================================================

body_raw <- read.csv(file.path(fit_dir, "body_stats.csv"),
                     header = TRUE, stringsAsFactors = FALSE)
cat("Body stats loaded:", nrow(body_raw), "individuals\n")

# Identify non-numeric weight entries before coercion
non_numeric_weight <- body_raw$ID[suppressWarnings(is.na(as.numeric(body_raw$weight))) &
                                    !is.na(body_raw$weight) & body_raw$weight != ""]
if (length(non_numeric_weight) > 0) {
  cat("\nNon-numeric weight values (will be excluded):\n")
  print(body_raw[body_raw$ID %in% non_numeric_weight,
                 c("ID", "weight", "svl", "site")])
}

body <- body_raw %>%
  mutate(
    ID     = as.character(ID),
    weight = suppressWarnings(as.numeric(weight)),
    svl    = as.numeric(svl),
    smi    = as.numeric(smi),
    site   = as.factor(site)
  ) %>%
  filter(!is.na(smi))   # require valid pre-calculated SMI

cat("After excluding missing SMI:", nrow(body), "individuals\n")

# =============================================================================
# 3. DERIVE BOUNDED SMI FROM PRE-CALCULATED SMI
# Scales pre-calculated SMI to (0,1) using Smithson & Verkuilen (2006)
# compression so Beta regression can be applied
# =============================================================================

n_obs   <- nrow(body)
smi_min <- min(body$smi, na.rm = TRUE)
smi_max <- max(body$smi, na.rm = TRUE)

body <- body %>%
  mutate(
    smi_scaled  = (smi - smi_min) / (smi_max - smi_min),
    smi_bounded = (smi_scaled * (n_obs - 1) + 0.5) / n_obs
  )

cat("SMI range:", round(smi_min, 3), "to", round(smi_max, 3), "\n")
cat("Bounded SMI range:", round(min(body$smi_bounded), 4),
    "to", round(max(body$smi_bounded), 4), "\n")

# =============================================================================
# 4. MERGE
# =============================================================================

hetfit <- inner_join(het_clean, body, by = "ID") %>%
  mutate(het = heterozygosity)

cat("\nMerged dataset:", nrow(hetfit), "individuals\n")
cat("By site:\n")
print(table(hetfit$site))

# =============================================================================
# 5. HANDLE SINGLETON SITES
# Daniel's Head (n=1) creates perfect leverage (hat=1) in site models,
# making HC3 robust SEs undefined and inflating interaction terms.
# Excluded from site models; retained in het-only and null models.
# =============================================================================

site_counts <- table(hetfit$site)
singleton_sites <- names(site_counts[site_counts < 3])
cat("\nSingleton/very small sites excluded from site models:", singleton_sites, "\n")

hetfit_site <- hetfit %>%
  filter(!site %in% singleton_sites) %>%
  mutate(site = droplevels(site))

cat("n for site models:", nrow(hetfit_site), "\n")

# =============================================================================
# 6. DIAGNOSTICS
# =============================================================================

cat("\n=== Normality checks ===\n")
cat("Shapiro-Wilk SMI:           p =", round(shapiro.test(hetfit$smi)$p.value, 4), "\n")
cat("Shapiro-Wilk log(SMI):      p =", round(shapiro.test(log(hetfit$smi))$p.value, 4), "\n")
cat("Shapiro-Wilk smi_bounded:   p =", round(shapiro.test(hetfit$smi_bounded)$p.value, 4), "\n")
cat("Shapiro-Wilk het:           p =", round(shapiro.test(hetfit$het)$p.value, 4), "\n")

cat("\n=== Summary by site ===\n")
hetfit %>%
  group_by(site) %>%
  summarise(n        = n(),
            mean_het = round(mean(het), 6),
            sd_het   = round(sd(het), 6),
            mean_smi = round(mean(smi, na.rm = TRUE), 4),
            sd_smi   = round(sd(smi, na.rm = TRUE), 4),
            .groups  = "drop") %>%
  print()

# =============================================================================
# 7. APPROACH 1 — OLS ON SMI WITH ROBUST STANDARD ERRORS (HC1)
# HC1 = (n/(n-k)) * HC0; stable when hat values approach 1
# Daniel's Head excluded from site models (singleton site)
# =============================================================================

cat("\n\n=== APPROACH 1: OLS on SMI — robust SEs (HC1) ===\n")

# Site models (singleton sites excluded)
r_full <- lm(smi ~ het * site, data = hetfit_site)
r_add  <- lm(smi ~ het + site, data = hetfit_site)

# Het-only and null use full dataset
r_het  <- lm(smi ~ het,        data = hetfit)
r_null <- lm(smi ~ 1,          data = hetfit)

cat("\nM_add (smi ~ het + site) — robust SEs (HC1):\n")
print(coeftest(r_add, vcov = vcovHC(r_add, type = "HC1")))

cat("\nWald test: interaction (M_full vs M_add):\n")
print(waldtest(r_full, r_add, vcov = vcovHC(r_full, type = "HC1")))

cat("\nWald test: site effect (M_add vs M_het):\n")
# Use F-test on M_add vs het-only fitted on same subset
r_het_sub <- lm(smi ~ het, data = hetfit_site)
print(waldtest(r_add, r_het_sub, vcov = vcovHC(r_add, type = "HC1")))

cat("\nWald test: het effect (M_het vs null — full dataset):\n")
print(waldtest(r_het, r_null, vcov = vcovHC(r_het, type = "HC1")))

cat("\nResidual normality M_add:", round(shapiro.test(residuals(r_add))$p.value, 4), "\n")

# =============================================================================
# 8. APPROACH 2 — OLS ON log(SMI)
# =============================================================================

cat("\n\n=== APPROACH 2: OLS on log(SMI) ===\n")

l_full <- lm(log(smi) ~ het * site, data = hetfit_site)
l_add  <- lm(log(smi) ~ het + site, data = hetfit_site)
l_het  <- lm(log(smi) ~ het,        data = hetfit)
l_null <- lm(log(smi) ~ 1,          data = hetfit)
l_het_sub <- lm(log(smi) ~ het,     data = hetfit_site)

cat("\nM_add log(smi) ~ het + site:\n")
print(summary(l_add))

cat("\nLRT interaction:\n");       print(anova(l_full, l_add))
cat("\nLRT site effect:\n");       print(anova(l_add, l_het_sub))
cat("\nLRT het effect (full n):\n"); print(anova(l_het, l_null))

cat("\nResidual normality log(SMI) M_add:",
    round(shapiro.test(residuals(l_add))$p.value, 4), "\n")

# =============================================================================
# 9. APPROACH 3 — BETA REGRESSION ON BOUNDED SMI
# Full interaction model dropped — overparameterised with small n per site
# =============================================================================

cat("\n\n=== APPROACH 3: Beta regression on bounded SMI ===\n")

hf_beta      <- hetfit      %>% filter(smi_bounded > 0, smi_bounded < 1)
hf_beta_site <- hetfit_site %>% filter(smi_bounded > 0, smi_bounded < 1)

cat("n for Beta additive model (site):", nrow(hf_beta_site), "\n")
cat("n for Beta het-only model:", nrow(hf_beta), "\n")

b_add  <- betareg(smi_bounded ~ het + site | 1, data = hf_beta_site)
b_het  <- betareg(smi_bounded ~ het        | 1, data = hf_beta)
b_null <- betareg(smi_bounded ~ 1          | 1, data = hf_beta)
b_het_sub <- betareg(smi_bounded ~ het     | 1, data = hf_beta_site)

cat("\nBeta M_add (het + site):\n")
print(summary(b_add))

cat("\nLRT site effect (Beta):\n")
print(lrtest(b_add, b_het_sub))

cat("\nLRT het effect (Beta, full n):\n")
print(lrtest(b_het, b_null))

cat("\nPseudo-R2 Beta M_add:", round(b_add$pseudo.r.squared, 4), "\n")

# =============================================================================
# 10. SUMMARY COMPARISON TABLE
# =============================================================================

cat("\n\n=== MODEL COMPARISON SUMMARY ===\n")

safe_p <- function(x) tryCatch(round(x, 4), error = function(e) NA)

comparison <- data.frame(
  approach      = c("OLS robust HC1 (SMI)",
                    "OLS (log SMI)",
                    "Beta regression"),
  n_site_model  = c(nrow(hetfit_site), nrow(hetfit_site), nrow(hf_beta_site)),
  het_coef      = round(c(coef(r_add)["het"],
                          coef(l_add)["het"],
                          coef(b_add)["het"]), 4),
  het_p_site_model = safe_p(c(
    coeftest(r_add, vcov = vcovHC(r_add, "HC1"))["het", "Pr(>|t|)"],
    summary(l_add)$coefficients["het", "Pr(>|t|)"],
    summary(b_add)$coefficients$mean["het", "Pr(>|z|)"]
  )),
  het_p_het_only = safe_p(c(
    waldtest(r_het, r_null, vcov = vcovHC(r_het, "HC1"))$`Pr(>F)`[2],
    anova(l_het, l_null)$`Pr(>F)`[2],
    lrtest(b_het, b_null)$`Pr(>Chisq)`[2]
  )),
  site_p = safe_p(c(
    waldtest(r_add, r_het_sub, vcov = vcovHC(r_add, "HC1"))$`Pr(>F)`[2],
    anova(l_add, l_het_sub)$`Pr(>F)`[2],
    lrtest(b_add, b_het_sub)$`Pr(>Chisq)`[2]
  )),
  interaction_p = safe_p(c(
    waldtest(r_full, r_add, vcov = vcovHC(r_full, "HC1"))$`Pr(>F)`[2],
    anova(l_full, l_add)$`Pr(>F)`[2],
    NA  # overparameterised for Beta
  ))
)

print(comparison)
write.table(comparison, file.path(out_dir, "hetfit_model_comparison.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)
cat("Saved: hetfit_model_comparison.tsv\n")

# =============================================================================
# 11. PLOTS
# =============================================================================

p_het <- ggplot(hetfit, aes(x = site, y = het, fill = site, colour = site)) +
  geom_violin(alpha = 0.2, linewidth = 0.4) +
  geom_boxplot(width = 0.12, outlier.shape = NA, fill = "white") +
  geom_jitter(width = 0.1, size = 1.5, alpha = 0.6) +
  scale_fill_manual(values = pop_colours, guide = "none") +
  scale_colour_manual(values = pop_colours, guide = "none") +
  labs(title = "Individual heterozygosity by site",
       x = NULL, y = "Per-site heterozygosity") +
  theme_bw(base_size = 11) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1),
        panel.grid.minor = element_blank())

p_smi <- ggplot(hetfit, aes(x = het, y = smi, colour = site, fill = site)) +
  geom_point(size = 2, alpha = 0.7) +
  geom_smooth(data = hetfit_site,
              method = "lm", se = TRUE, alpha = 0.12, linewidth = 0.8) +
  scale_colour_manual(values = pop_colours, name = "Site") +
  scale_fill_manual(values = pop_colours, name = "Site") +
  labs(title = "Heterozygosity vs SMI by site",
       x = "Per-site heterozygosity", y = "SMI (g)") +
  theme_bw(base_size = 11) +
  theme(panel.grid.minor = element_blank(),
        legend.position = "bottom",
        legend.text = element_text(size = 8))

p_logsmi <- ggplot(hetfit, aes(x = het, y = log(smi))) +
  geom_point(aes(colour = site), size = 2, alpha = 0.7) +
  geom_smooth(method = "lm", se = TRUE,
              colour = "grey20", fill = "grey70", alpha = 0.2) +
  scale_colour_manual(values = pop_colours, name = "Site") +
  labs(title = "Heterozygosity vs log(SMI) — overall trend",
       x = "Per-site heterozygosity", y = "log(SMI)") +
  theme_bw(base_size = 11) +
  theme(panel.grid.minor = element_blank(),
        legend.position = "bottom",
        legend.text = element_text(size = 8))

p_beta <- ggplot(hf_beta, aes(x = het, y = smi_bounded)) +
  geom_point(aes(colour = site), size = 2, alpha = 0.7) +
  geom_smooth(method = "lm", se = TRUE,
              colour = "grey20", fill = "grey70", alpha = 0.2) +
  scale_colour_manual(values = pop_colours, name = "Site") +
  labs(title = "Heterozygosity vs bounded SMI (Beta response)",
       x = "Per-site heterozygosity", y = "Bounded SMI (0–1)") +
  theme_bw(base_size = 11) +
  theme(panel.grid.minor = element_blank(),
        legend.position = "bottom",
        legend.text = element_text(size = 8))

qq_data <- bind_rows(
  data.frame(residuals = residuals(r_add), model = "OLS: SMI ~ het + site"),
  data.frame(residuals = residuals(l_add), model = "OLS: log(SMI) ~ het + site"),
  data.frame(residuals = residuals(b_add), model = "Beta: bounded SMI ~ het + site")
)

p_qq <- ggplot(qq_data, aes(sample = residuals)) +
  stat_qq(size = 0.8, alpha = 0.6) +
  stat_qq_line(colour = "firebrick", linewidth = 0.6) +
  facet_wrap(~model, scales = "free") +
  labs(title = "Residual Q-Q plots — additive models",
       x = "Theoretical quantiles", y = "Sample quantiles") +
  theme_bw(base_size = 11) +
  theme(panel.grid.minor = element_blank())

p_combined <- (p_het | p_logsmi) /
  (p_smi | p_beta) /
  p_qq +
  plot_annotation(
    title    = "Heterozygosity-fitness correlations — Plestiodon longirostris",
    subtitle = paste0("Three approaches: robust OLS (HC1), log(SMI) OLS, Beta regression  |  ",
                      "Daniel's Head (n=1) excluded from site models"),
    theme    = theme(plot.title    = element_text(size = 13, face = "bold"),
                     plot.subtitle = element_text(size = 9, colour = "grey40"))
  ) +
  plot_layout(heights = c(1, 1, 1))

ggsave(file.path(out_dir, "hetfit_combined.pdf"),
       p_combined, width = 14, height = 16)
ggsave(file.path(out_dir, "hetfit_combined.png"),
       p_combined, width = 14, height = 16, dpi = 300)
cat("Saved: hetfit_combined.pdf/.png\n")

ggsave(file.path(out_dir, "hetfit_het_dist.pdf"),  p_het,    width = 8, height = 5)
ggsave(file.path(out_dir, "hetfit_smi_site.pdf"),  p_smi,    width = 9, height = 6)
ggsave(file.path(out_dir, "hetfit_logsmi.pdf"),    p_logsmi, width = 7, height = 6)
ggsave(file.path(out_dir, "hetfit_beta.pdf"),      p_beta,   width = 7, height = 6)
ggsave(file.path(out_dir, "hetfit_qq.pdf"),        p_qq,     width = 11, height = 4)

# =============================================================================
# 12. POST-HOC SITE COMPARISONS (emmeans)
# Pairwise Tukey-corrected comparisons of site marginal means
# Based on log(SMI) additive model — controls for heterozygosity
# emmeans gives estimated marginal means at the mean heterozygosity
# =============================================================================

if (!requireNamespace("emmeans",  quietly = TRUE)) install.packages("emmeans")
if (!requireNamespace("multcomp", quietly = TRUE)) install.packages("multcomp")
library(emmeans)
library(multcomp)

cat("\n\n=== POST-HOC SITE COMPARISONS — log(SMI) additive model ===\n")
cat("Marginal means estimated at mean heterozygosity\n\n")

emm <- emmeans(l_add, ~ site)

cat("Estimated marginal means per site:\n")
print(emm)

cat("\nPairwise contrasts (Tukey correction):\n")
emm_pairs <- pairs(emm, adjust = "tukey")
print(emm_pairs)

# Save pairwise table
emm_df <- as.data.frame(emm_pairs)
write.table(emm_df,
            file.path(out_dir, "hetfit_site_pairwise_tukey.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)
cat("Saved: hetfit_site_pairwise_tukey.tsv\n")

# Compact letter display via cld() is version-dependent and unreliable.
# Use pairwise significance table instead — equivalent information.
cat("\nPairwise significance (Tukey corrected):\n")
emm_sig <- as.data.frame(emm_pairs) %>%
  mutate(significant = ifelse(p.value < 0.05, "*", "ns")) %>%
  dplyr::select(contrast, estimate, SE, t.ratio, p.value, significant)
print(emm_sig)

cld_df <- as.data.frame(emm) %>%
  dplyr::select(site, emmean, lower.CL, upper.CL) %>%
  mutate(site = as.factor(site))

# --- Emmeans plot: site marginal means with comparison arrows ---
p_emm <- plot(emm, comparisons = TRUE) +
  labs(
    title    = "Site marginal means — log(SMI)",
    subtitle = paste0("Estimated at mean heterozygosity  |  ",
                      "Arrows: non-overlapping = significantly different (Tukey)"),
    x        = "log(SMI) marginal mean",
    y        = "Site"
  ) +
  theme_bw(base_size = 12) +
  theme(panel.grid.minor = element_blank())

ggsave(file.path(out_dir, "hetfit_site_emmeans.pdf"),
       p_emm, width = 8, height = 5)
ggsave(file.path(out_dir, "hetfit_site_emmeans.png"),
       p_emm, width = 8, height = 5, dpi = 300)
cat("Saved: hetfit_site_emmeans.pdf/.png\n")

# --- SMI violin by site with significant pairs annotated ---
# Mark sites that differ significantly from Castle Island (reference)
sig_vs_ref <- emm_sig %>%
  mutate(
    site = ifelse(grepl("^Castle Island", contrast),
                  sub("Castle Island - ", "", contrast),
                  ifelse(grepl("Castle Island$", contrast),
                         sub(" - Castle Island", "", contrast),
                         NA_character_)),
    label = ifelse(significant == "*" & !is.na(site), "*", "")
  ) %>%
  filter(!is.na(site))

site_max <- hetfit %>%
  group_by(site) %>%
  summarise(y_pos = max(smi, na.rm = TRUE) + 0.025, .groups = "drop") %>%
  mutate(site = as.character(site)) %>%
  left_join(sig_vs_ref %>% dplyr::select(site, label), by = "site") %>%
  mutate(label = ifelse(is.na(label), "", label),
         label = ifelse(site == "Castle Island", "ref", label))

p_smi_sig <- ggplot(hetfit,
                    aes(x = site, y = smi, fill = site, colour = site)) +
  geom_violin(alpha = 0.2, linewidth = 0.4) +
  geom_boxplot(width = 0.12, outlier.shape = NA,
               fill = "white", linewidth = 0.5) +
  geom_jitter(width = 0.1, size = 1.2, alpha = 0.5) +
  geom_text(data = site_max,
            aes(x = site, y = y_pos, label = label),
            inherit.aes = FALSE,
            size = 4.5, fontface = "bold", colour = "firebrick") +
  scale_fill_manual(values = pop_colours, guide = "none") +
  scale_colour_manual(values = pop_colours, guide = "none") +
  labs(
    title    = "SMI by site — post-hoc comparisons",
    subtitle = "* = significantly different from Castle Island (Tukey p < 0.05)  |  ref = reference site",
    x        = NULL,
    y        = "SMI (g)"
  ) +
  theme_bw(base_size = 12) +
  theme(axis.text.x    = element_text(angle = 30, hjust = 1),
        panel.grid.minor = element_blank())

ggsave(file.path(out_dir, "hetfit_smi_significance.pdf"),
       p_smi_sig, width = 9, height = 6)
ggsave(file.path(out_dir, "hetfit_smi_significance.png"),
       p_smi_sig, width = 9, height = 6, dpi = 300)
cat("Saved: hetfit_smi_significance.pdf/.png\n")

# =============================================================================
# 13. SAVE MERGED DATA
# =============================================================================

write.table(
  hetfit %>%
    select(ID, site, het, n_sites_total, weight, svl,
           smi, smi_bounded) %>%
    arrange(site, ID),
  file.path(out_dir, "hetfit_merged_data.tsv"),
  sep = "\t", row.names = FALSE, quote = FALSE
)

cat("\nAll outputs saved to:", out_dir, "\n")
