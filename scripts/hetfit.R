# =============================================================================
# Heterozygosity-Fitness Regression
# Plestiodon longirostris - Bermuda lcWGS
#
# SMI calculated from scratch using bSMA = SD(ln mass) / SD(ln SVL)
# Sites with n < 5 excluded from all models
#
# Model selection:
#   Three candidate models fitted (OLS, log OLS, Beta regression)
#   Best model selected automatically by:
#     1. Residual normality (Shapiro-Wilk p > 0.05 preferred)
#     2. Outlier influence (Cook's D — flag if any > 4/n)
#     3. AIC (lower = better fit)
#   Best model results reported; others summarised for comparison
#
# Reference: Peig & Green (2009) Oikos 118:1883-1891
# =============================================================================

library(dplyr)
library(tidyr)
library(ggplot2)
library(patchwork)
library(sandwich)
library(lmtest)
library(betareg)

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

MIN_SITE_N <- 5   # minimum n per site for inclusion in site models
MIN_SVL    <- 65  # minimum SVL (mm) — exclude subadults from analysis
# Distribution shows clear break at ~65mm between subadults
# and adults (main peak ~80mm). Subadults have different
# bSMA slope (growth-related mass gain) vs adults (condition).

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
# 2. LOAD BODY STATS AND CALCULATE SMI FROM SCRATCH
# bSMA = SD(ln mass) / SD(ln SVL)  [standardised major axis slope]
# SMI  = mass * (mean_SVL / SVL) ^ bSMA  [Peig & Green 2009]
# =============================================================================

body_raw <- read.csv(file.path(fit_dir, "body_stats.csv"),
                     header = TRUE, stringsAsFactors = FALSE)
cat("Body stats loaded:", nrow(body_raw), "individuals\n")

non_numeric_weight <- body_raw$ID[
  suppressWarnings(is.na(as.numeric(body_raw$weight))) &
    !is.na(body_raw$weight) & body_raw$weight != ""]
if (length(non_numeric_weight) > 0) {
  cat("\nNon-numeric weight values excluded from bSMA calculation:\n")
  print(body_raw[body_raw$ID %in% non_numeric_weight,
                 c("ID", "weight", "svl", "site")])
}

body <- body_raw %>%
  mutate(
    ID     = as.character(ID),
    weight = suppressWarnings(as.numeric(weight)),
    svl    = as.numeric(svl),
    site   = as.factor(site)
  )

body_valid <- body %>%
  filter(!is.na(weight), !is.na(svl), weight > 0, svl >= MIN_SVL) %>%
  mutate(ln_mass = log(weight), ln_svl = log(svl))

bSMA <- sd(body_valid$ln_mass) / sd(body_valid$ln_svl)
L0   <- mean(body_valid$svl)

cat("\n=== SMI parameters (adults only, SVL >=", MIN_SVL, "mm) ===\n")
cat("bSMA:", round(bSMA, 4), " | L0:", round(L0, 2), "mm |",
    nrow(body_valid), "individuals used\n")

if (bSMA < 2 || bSMA > 4) {
  cat("WARNING: bSMA outside expected range (2-4) for lizards\n")
} else {
  cat("bSMA within expected range for lizards (2-4) — OK\n")
}

body <- body %>%
  mutate(smi = ifelse(!is.na(weight) & !is.na(svl) & weight > 0 & svl >= MIN_SVL,
                      weight * (L0 / svl) ^ bSMA, NA_real_)) %>%
  filter(!is.na(smi))

cat("Individuals with valid SMI (SVL >=", MIN_SVL, "mm):", nrow(body), "\n")
cat("SMI range:", round(min(body$smi), 2), "to", round(max(body$smi), 2),
    "g | Mean:", round(mean(body$smi), 2), "g\n")

# Bounded SMI for Beta regression
n_obs <- nrow(body)
body <- body %>%
  mutate(
    smi_scaled  = (smi - min(smi)) / (max(smi) - min(smi)),
    smi_bounded = (smi_scaled * (n_obs - 1) + 0.5) / n_obs
  )

# =============================================================================
# 3. MERGE AND FILTER SITES
# Exclude sites with n < MIN_SITE_N from all analyses
# =============================================================================

hetfit_all <- inner_join(het_clean, body, by = "ID") %>%
  mutate(het = heterozygosity)

cat("\nAll merged individuals:", nrow(hetfit_all), "\n")
cat("By site (before filtering):\n")
print(table(hetfit_all$site))

site_counts <- table(hetfit_all$site)
small_sites <- names(site_counts[site_counts < MIN_SITE_N])
cat("\nSites excluded (n <", MIN_SITE_N, "):", paste(small_sites, collapse = ", "), "\n")

hetfit <- hetfit_all %>%
  filter(!site %in% small_sites) %>%
  mutate(site = droplevels(site))

hf_beta <- hetfit %>% filter(smi_bounded > 0, smi_bounded < 1)

cat("Individuals retained:", nrow(hetfit), "\n")
cat("By site:\n")
print(table(hetfit$site))

# =============================================================================
# 4. DIAGNOSTICS
# =============================================================================

cat("\n=== Summary by site ===\n")
hetfit %>%
  group_by(site) %>%
  summarise(n        = n(),
            mean_het = round(mean(het), 6),
            sd_het   = round(sd(het), 6),
            mean_smi = round(mean(smi), 3),
            sd_smi   = round(sd(smi), 3),
            .groups  = "drop") %>%
  print()

# =============================================================================
# 5. FIT ALL THREE CANDIDATE MODELS
# Each model: additive (het + site), het-only, null
# =============================================================================

cat("\n=== Fitting candidate models ===\n")

# OLS on SMI
r_add  <- lm(smi       ~ het + site, data = hetfit)
r_full <- lm(smi       ~ het * site, data = hetfit)
r_het  <- lm(smi       ~ het,        data = hetfit)
r_null <- lm(smi       ~ 1,          data = hetfit)

# OLS on log(SMI)
l_add  <- lm(log(smi)  ~ het + site, data = hetfit)
l_full <- lm(log(smi)  ~ het * site, data = hetfit)
l_het  <- lm(log(smi)  ~ het,        data = hetfit)
l_null <- lm(log(smi)  ~ 1,          data = hetfit)

# OLS on sqrt(SMI)
s_add  <- lm(sqrt(smi) ~ het + site, data = hetfit)
s_full <- lm(sqrt(smi) ~ het * site, data = hetfit)
s_het  <- lm(sqrt(smi) ~ het,        data = hetfit)
s_null <- lm(sqrt(smi) ~ 1,          data = hetfit)

# Beta regression
b_add     <- betareg(smi_bounded ~ het + site | 1, data = hf_beta)
b_het     <- betareg(smi_bounded ~ het        | 1, data = hf_beta)
b_null    <- betareg(smi_bounded ~ 1          | 1, data = hf_beta)

# =============================================================================
# 6. MODEL DIAGNOSTICS — select best model automatically
#
# Criteria (applied to additive model of each type):
#   A. Residual normality: Shapiro-Wilk p-value (higher = better)
#   B. Outlier influence: max Cook's D relative to threshold 4/n
#      (for Beta, use standardised residuals instead)
#   C. AIC (lower = better; Beta AIC from log-likelihood)
#
# Scoring: each model gets 0-3 points
#   +1 if SW p > 0.05 (residuals normal)
#   +1 if no influential outliers (max CookD < 4/n)
#   +1 if lowest AIC of the three
# Model with most points selected; ties broken by AIC
# =============================================================================

cat("\n=== Model diagnostics ===\n")

n           <- nrow(hetfit)
n_beta      <- nrow(hf_beta)
cook_thresh <- 4 / n

# Residual normality (Shapiro-Wilk)
sw_ols  <- shapiro.test(residuals(r_add))$p.value
sw_log  <- shapiro.test(residuals(l_add))$p.value
sw_sqrt <- shapiro.test(residuals(s_add))$p.value
sw_beta <- shapiro.test(residuals(b_add, type = "quantile"))$p.value

# Cook's distance for OLS models
max_cook_ols  <- max(cooks.distance(r_add), na.rm = TRUE)
max_cook_log  <- max(cooks.distance(l_add), na.rm = TRUE)
max_cook_sqrt <- max(cooks.distance(s_add), na.rm = TRUE)

# Beta: quantile residual outliers (|z| > 3)
q_beta     <- residuals(b_add, type = "quantile")
outlier_beta <- max(abs(q_beta), na.rm = TRUE) > 3.0

# AIC
aic_ols  <- AIC(r_add)
aic_log  <- AIC(l_add)
aic_sqrt <- AIC(s_add)
aic_beta <- as.numeric(-2 * logLik(b_add) + 2 * (length(coef(b_add)) + 1))

cat(sprintf("%-26s %8s %12s %10s %10s\n",
            "Model", "SW p", "Max Cook's D", "AIC", "Outliers?"))
cat(sprintf("%-26s %8.4f %12.4f %10.2f %10s\n",
            "OLS (SMI)", sw_ols, max_cook_ols, aic_ols,
            ifelse(max_cook_ols  > cook_thresh, "YES", "no")))
cat(sprintf("%-26s %8.4f %12.4f %10.2f %10s\n",
            "OLS (log SMI)", sw_log, max_cook_log, aic_log,
            ifelse(max_cook_log  > cook_thresh, "YES", "no")))
cat(sprintf("%-26s %8.4f %12.4f %10.2f %10s\n",
            "OLS (sqrt SMI)", sw_sqrt, max_cook_sqrt, aic_sqrt,
            ifelse(max_cook_sqrt > cook_thresh, "YES", "no")))
cat(sprintf("%-26s %8.4f %12s %10.2f %10s\n",
            "Beta regression", sw_beta, "(quantile)", aic_beta,
            ifelse(outlier_beta, "YES", "no")))

# Score: +1 normal residuals, +1 no outliers, +1 lowest AIC
score_ols  <- (sw_ols  > 0.05) + (max_cook_ols  <= cook_thresh)
score_log  <- (sw_log  > 0.05) + (max_cook_log  <= cook_thresh)
score_sqrt <- (sw_sqrt > 0.05) + (max_cook_sqrt <= cook_thresh)
score_beta <- (sw_beta > 0.05) + (!outlier_beta)

aics <- c(OLS = aic_ols, Log_OLS = aic_log,
          Sqrt_OLS = aic_sqrt, Beta = aic_beta)
score_ols  <- score_ols  + (which.min(aics) == 1)
score_log  <- score_log  + (which.min(aics) == 2)
score_sqrt <- score_sqrt + (which.min(aics) == 3)
score_beta <- score_beta + (which.min(aics) == 4)

scores <- c(OLS = score_ols, Log_OLS = score_log,
            Sqrt_OLS = score_sqrt, Beta = score_beta)
cat("\nModel scores (max 3):",
    paste(names(scores), scores, sep = "=", collapse = "  "), "\n")

# Strip any "name.name" prefix R sometimes adds when subsetting named vectors
clean_name <- function(x) sub("^.*\\.", "", x)

# Best = highest score; ties broken by best SW p (most normal residuals)
best_name <- clean_name(names(scores)[which.max(scores)])
if (sum(scores == max(scores)) > 1) {
  tied    <- clean_name(names(scores[scores == max(scores)]))
  sw_vals <- setNames(c(sw_ols, sw_log, sw_sqrt, sw_beta),
                      c("OLS", "Log_OLS", "Sqrt_OLS", "Beta"))[tied]
  best_name <- names(sw_vals)[which.max(sw_vals)]
  cat("Tie broken by best SW p — best model:", best_name, "\n")
}
cat("==> Best model selected:", best_name, "\n\n")
# =============================================================================
# 7. REPORT BEST MODEL FULLY
# =============================================================================

run_ols_report <- function(add_mod, full_mod, het_mod, null_mod, response_label) {
  cat("--- Additive model:", response_label, "~ het + site ---\n")
  print(summary(add_mod))
  
  cat("\nLRT: interaction (", response_label, "~ het*site vs het+site):\n")
  print(anova(full_mod, add_mod))
  
  cat("\nLRT: site effect (het+site vs het):\n")
  print(anova(add_mod, het_mod))
  
  cat("\nLRT: het effect (het vs null):\n")
  print(anova(het_mod, null_mod))
  
  cat("\nResidual normality:", round(shapiro.test(residuals(add_mod))$p.value, 4), "\n")
  
  n_cook <- length(cooks.distance(add_mod))
  thresh  <- 4 / n_cook
  n_inf   <- sum(cooks.distance(add_mod) > thresh)
  cat("Influential observations (Cook's D >", round(thresh, 4), "):", n_inf, "\n")
  if (n_inf > 0) {
    inf_ids <- hetfit$ID[cooks.distance(add_mod) > thresh]
    cat("  IDs:", paste(inf_ids, collapse = ", "), "\n")
  }
}

run_robust_report <- function(add_mod, full_mod, het_sub, null_mod, n) {
  cat("--- Additive model: SMI ~ het + site (HC1 robust SEs) ---\n")
  print(coeftest(add_mod, vcov = vcovHC(add_mod, "HC1")))
  
  cat("\nWald test: interaction:\n")
  print(waldtest(full_mod, add_mod, vcov = vcovHC(full_mod, "HC1")))
  
  cat("\nWald test: site effect:\n")
  print(waldtest(add_mod, het_sub, vcov = vcovHC(add_mod, "HC1")))
  
  cat("\nWald test: het effect:\n")
  het_only <- lm(smi ~ het, data = hetfit)
  null_mod2 <- lm(smi ~ 1,  data = hetfit)
  print(waldtest(het_only, null_mod2, vcov = vcovHC(het_only, "HC1")))
  
  cat("\nResidual normality:", round(shapiro.test(residuals(add_mod))$p.value, 4), "\n")
  thresh <- 4 / n
  n_inf  <- sum(cooks.distance(add_mod) > thresh)
  cat("Influential observations (Cook's D >", round(thresh, 4), "):", n_inf, "\n")
  if (n_inf > 0) {
    inf_ids <- hetfit$ID[cooks.distance(add_mod) > thresh]
    cat("  IDs:", paste(inf_ids, collapse = ", "), "\n")
  }
}

run_beta_report <- function(add_mod, het_mod, het_sub, null_mod) {
  cat("--- Beta regression: bounded SMI ~ het + site ---\n")
  print(summary(add_mod))
  
  cat("\nLRT: site effect:\n")
  print(lrtest(add_mod, het_sub))
  
  cat("\nLRT: het effect:\n")
  print(lrtest(het_mod, null_mod))
  
  cat("\nPseudo-R2:", round(add_mod$pseudo.r.squared, 4), "\n")
  
  q_resid <- residuals(add_mod, type = "quantile")
  cat("Residual normality (quantile):",
      round(shapiro.test(q_resid)$p.value, 4), "\n")
  cat("Outliers (|quantile resid| > 3):",
      sum(abs(q_resid) > 3), "\n")
}

cat("=============================================================\n")
cat("BEST MODEL RESULTS:", best_name, "\n")
cat("=============================================================\n\n")

r_het_sub <- lm(smi       ~ het, data = hetfit)
l_het_sub <- lm(log(smi)  ~ het, data = hetfit)
s_het_sub <- lm(sqrt(smi) ~ het, data = hetfit)
b_het_sub <- betareg(smi_bounded ~ het | 1, data = hf_beta)

if (best_name == "OLS") {
  run_robust_report(r_add, r_full, r_het_sub, r_null, n)
} else if (best_name == "Log_OLS") {
  run_ols_report(l_add, l_full, l_het_sub, l_null, "log(SMI)")
} else if (best_name == "Sqrt_OLS") {
  run_ols_report(s_add, s_full, s_het_sub, s_null, "sqrt(SMI)")
} else {
  run_beta_report(b_add, b_het, b_het_sub, b_null)
}

# =============================================================================
# 8. COMPARISON SUMMARY TABLE (all three models)
# =============================================================================

cat("\n\n=== FULL MODEL COMPARISON (all three approaches) ===\n")

safe_p <- function(x) tryCatch(round(as.numeric(x), 4), error = function(e) NA)

comparison <- data.frame(
  approach      = c("OLS HC1 (SMI)", "OLS log(SMI)", "OLS sqrt(SMI)", "Beta regression"),
  selected      = c(best_name == "OLS", best_name == "Log_OLS",
                    best_name == "Sqrt_OLS", best_name == "Beta"),
  n             = c(nrow(hetfit), nrow(hetfit), nrow(hetfit), nrow(hf_beta)),
  SW_p          = round(c(sw_ols, sw_log, sw_sqrt, sw_beta), 4),
  max_CooksD    = round(c(max_cook_ols, max_cook_log, max_cook_sqrt, NA), 4),
  AIC           = round(c(aic_ols, aic_log, aic_sqrt, aic_beta), 2),
  score         = c(score_ols, score_log, score_sqrt, score_beta),
  het_p_additive = safe_p(c(
    coeftest(r_add, vcov = vcovHC(r_add, "HC1"))["het", "Pr(>|t|)"],
    summary(l_add)$coefficients["het", "Pr(>|t|)"],
    summary(s_add)$coefficients["het", "Pr(>|t|)"],
    summary(b_add)$coefficients$mean["het", "Pr(>|z|)"]
  )),
  het_p_het_only = safe_p(c(
    waldtest(r_het_sub, lm(smi ~ 1, data = hetfit),
             vcov = vcovHC(r_het_sub, "HC1"))$`Pr(>F)`[2],
    anova(l_het, l_null)$`Pr(>F)`[2],
    anova(s_het, s_null)$`Pr(>F)`[2],
    lrtest(b_het, b_null)$`Pr(>Chisq)`[2]
  )),
  site_p = safe_p(c(
    waldtest(r_add, r_het_sub, vcov = vcovHC(r_add, "HC1"))$`Pr(>F)`[2],
    anova(l_add, l_het_sub)$`Pr(>F)`[2],
    anova(s_add, s_het_sub)$`Pr(>F)`[2],
    lrtest(b_add, b_het_sub)$`Pr(>Chisq)`[2]
  )),
  interaction_p = safe_p(c(
    waldtest(r_full, r_add, vcov = vcovHC(r_full, "HC1"))$`Pr(>F)`[2],
    anova(l_full, l_add)$`Pr(>F)`[2],
    anova(s_full, s_add)$`Pr(>F)`[2],
    NA
  ))
)

print(comparison)
write.table(comparison,
            file.path(out_dir, "hetfit_model_comparison.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)
cat("Saved: hetfit_model_comparison.tsv\n")

# =============================================================================
# 9. POST-HOC SITE COMPARISONS (emmeans on best linear model)
# =============================================================================

if (!requireNamespace("emmeans", quietly = TRUE)) install.packages("emmeans")
library(emmeans)

best_lm <- if (best_name == "OLS")      r_add else
  if (best_name == "Log_OLS")  l_add else
    if (best_name == "Sqrt_OLS") s_add else
      l_add
best_label <- if (best_name == "OLS")      "SMI" else
  if (best_name == "Log_OLS")  "log(SMI)" else
    if (best_name == "Sqrt_OLS") "sqrt(SMI)" else
      "log(SMI)"

cat("\n\n=== POST-HOC SITE COMPARISONS ===\n")
cat("Model:", best_label, "~ het + site  |  Tukey correction\n\n")

emm       <- emmeans(best_lm, ~ site)
emm_pairs <- pairs(emm, adjust = "tukey")

cat("Estimated marginal means:\n")
print(emm)
cat("\nPairwise contrasts:\n")
print(emm_pairs)

emm_sig <- as.data.frame(emm_pairs) %>%
  mutate(significant = ifelse(p.value < 0.05, "*", "ns")) %>%
  select(contrast, estimate, SE, t.ratio, p.value, significant)

write.table(as.data.frame(emm_pairs),
            file.path(out_dir, "hetfit_site_pairwise_tukey.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)
cat("Saved: hetfit_site_pairwise_tukey.tsv\n")

cld_df <- as.data.frame(emm) %>%
  select(site, emmean, lower.CL, upper.CL) %>%
  mutate(site = as.factor(site))

# =============================================================================
# 10. PLOTS
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
  geom_smooth(method = "lm", se = TRUE, alpha = 0.12, linewidth = 0.8) +
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

# Diagnostic plots for best model
if (best_name %in% c("OLS", "Log_OLS", "Sqrt_OLS")) {
  best_mod <- if (best_name == "OLS") r_add else if (best_name == "Log_OLS") l_add else s_add
  diag_data <- data.frame(
    fitted    = fitted(best_mod),
    residuals = residuals(best_mod),
    cooks     = cooks.distance(best_mod),
    obs       = seq_len(nrow(hetfit))
  )
  
  p_resid <- ggplot(diag_data, aes(x = fitted, y = residuals)) +
    geom_point(alpha = 0.6, size = 1.5) +
    geom_hline(yintercept = 0, linetype = "dashed", colour = "firebrick") +
    geom_smooth(method = "loess", se = FALSE, colour = "steelblue",
                linewidth = 0.7) +
    labs(title = paste("Residuals vs fitted —", best_name),
         x = "Fitted values", y = "Residuals") +
    theme_bw(base_size = 11)
  
  p_qq_best <- ggplot(diag_data, aes(sample = residuals)) +
    stat_qq(size = 1, alpha = 0.7) +
    stat_qq_line(colour = "firebrick", linewidth = 0.6) +
    labs(title = paste("Q-Q plot —", best_name),
         x = "Theoretical", y = "Sample") +
    theme_bw(base_size = 11)
  
  p_cook <- ggplot(diag_data, aes(x = obs, y = cooks)) +
    geom_col(fill = "steelblue", alpha = 0.7) +
    geom_hline(yintercept = 4 / n, linetype = "dashed",
               colour = "firebrick", linewidth = 0.6) +
    annotate("text", x = max(diag_data$obs) * 0.95,
             y = 4 / n + max(diag_data$cooks) * 0.02,
             label = paste("threshold =", round(4/n, 3)),
             hjust = 1, size = 3, colour = "firebrick") +
    labs(title = paste("Cook's distance —", best_name),
         x = "Observation", y = "Cook's D") +
    theme_bw(base_size = 11)
  
  p_diag <- (p_resid | p_qq_best | p_cook) +
    plot_annotation(title = paste("Diagnostics: best model (", best_name, ")"))
} else {
  # Beta diagnostics
  q_resid <- residuals(b_add, type = "quantile")
  diag_beta <- data.frame(fitted = fitted(b_add), residuals = q_resid,
                          obs = seq_len(nrow(hf_beta)))
  p_resid <- ggplot(diag_beta, aes(x = fitted, y = residuals)) +
    geom_point(alpha = 0.6, size = 1.5) +
    geom_hline(yintercept = 0, linetype = "dashed", colour = "firebrick") +
    labs(title = "Residuals vs fitted — Beta", x = "Fitted", y = "Quantile residuals") +
    theme_bw(base_size = 11)
  p_qq_best <- ggplot(diag_beta, aes(sample = residuals)) +
    stat_qq(size = 1, alpha = 0.7) + stat_qq_line(colour = "firebrick") +
    labs(title = "Q-Q — Beta") + theme_bw(base_size = 11)
  p_diag <- (p_resid | p_qq_best) +
    plot_annotation(title = "Diagnostics: best model (Beta)")
}

# Emmeans plot
p_emm <- plot(emm, comparisons = TRUE) +
  labs(title = paste("Site marginal means —", best_label),
       subtitle = "Arrows: non-overlapping = significantly different (Tukey)",
       x = paste(best_label, "marginal mean"), y = "Site") +
  theme_bw(base_size = 12) +
  theme(panel.grid.minor = element_blank())

# SMI violin with significance annotations
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
  summarise(y_pos = max(smi, na.rm = TRUE) * 1.03, .groups = "drop") %>%
  mutate(site = as.character(site)) %>%
  left_join(sig_vs_ref %>% select(site, label), by = "site") %>%
  mutate(label = ifelse(is.na(label), "", label),
         label = ifelse(site == "Castle Island", "ref", label))

p_smi_sig <- ggplot(hetfit,
                    aes(x = site, y = smi, fill = site, colour = site)) +
  geom_violin(alpha = 0.2, linewidth = 0.4) +
  geom_boxplot(width = 0.12, outlier.shape = NA, fill = "white") +
  geom_jitter(width = 0.1, size = 1.2, alpha = 0.5) +
  geom_text(data = site_max,
            aes(x = site, y = y_pos, label = label),
            inherit.aes = FALSE, size = 5,
            fontface = "bold", colour = "firebrick") +
  scale_fill_manual(values = pop_colours, guide = "none") +
  scale_colour_manual(values = pop_colours, guide = "none") +
  labs(title = "SMI by site — post-hoc comparisons",
       subtitle = paste0("* = significantly different from Castle Island  |  ",
                         "Sites with n < ", MIN_SITE_N, " excluded"),
       x = NULL, y = "SMI (g)") +
  theme_bw(base_size = 12) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1),
        panel.grid.minor = element_blank())

# Combined panel
p_combined <- (p_het | p_logsmi) /
  (p_smi | p_emm) /
  p_diag +
  plot_annotation(
    title    = paste0("Heterozygosity-fitness correlations — Plestiodon longirostris",
                      "  |  Best model: ", best_name),
    subtitle = paste0("bSMA = ", round(bSMA, 3),
                      "  |  Sites with n < ", MIN_SITE_N, " excluded",
                      "  |  n = ", nrow(hetfit), " individuals"),
    theme    = theme(plot.title    = element_text(size = 13, face = "bold"),
                     plot.subtitle = element_text(size = 9, colour = "grey40"))
  ) +
  plot_layout(heights = c(1, 1, 1))

ggsave(file.path(out_dir, "hetfit_combined.pdf"),
       p_combined, width = 14, height = 16)
ggsave(file.path(out_dir, "hetfit_combined.png"),
       p_combined, width = 14, height = 16, dpi = 300)

ggsave(file.path(out_dir, "hetfit_smi_significance.pdf"), p_smi_sig, width = 9, height = 6)
ggsave(file.path(out_dir, "hetfit_smi_significance.png"), p_smi_sig, width = 9, height = 6, dpi = 300)
ggsave(file.path(out_dir, "hetfit_site_emmeans.pdf"),     p_emm,     width = 8, height = 5)
ggsave(file.path(out_dir, "hetfit_site_emmeans.png"),     p_emm,     width = 8, height = 5, dpi = 300)
ggsave(file.path(out_dir, "hetfit_smi_scatter.pdf"),      p_smi,     width = 9, height = 6)
ggsave(file.path(out_dir, "hetfit_het_dist.pdf"),         p_het,     width = 8, height = 5)
cat("Saved: all plots\n")

# =============================================================================
# 11. SAVE MERGED DATA
# =============================================================================

write.table(
  hetfit %>%
    select(ID, site, het, n_sites_total, weight, svl, smi, smi_bounded) %>%
    mutate(bSMA_used = round(bSMA, 6), L0_used = round(L0, 3)) %>%
    arrange(site, ID),
  file.path(out_dir, "hetfit_merged_data.tsv"),
  sep = "\t", row.names = FALSE, quote = FALSE
)

cat("\nAll outputs saved to:", out_dir, "\n")
cat("Best model:", best_name, "(score", max(scores), "/ 3)\n")

#For reporting you could state it as: "Heterozygosity was positively associated with scaled mass index 
#(β = 44.1, F₁,₁₅₄ = 4.94, p = 0.028), suggesting individuals with higher genome-wide heterozygosity were
#in better body condition, though this effect was partially attenuated when sampling site was included 
#as a covariate (p = 0.116)."