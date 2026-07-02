########################################################################
# Luttrell, Petty, & Brinol (2016) Study 1 — Other Analyses
# Includes:
#   Part A: Standardized regression (following rjb_analysis.R)
#   Part B: Three-way interaction  (following LuttrellInferentialTestsOSF.R)
#   Part C: Bayesian estimation with Bayes factor, HDI, & sensitivity
########################################################################

library(tidyverse)
library(psych)
library(Hmisc)
library(BayesFactor)
library(bayestestR)

select <- dplyr::select
library(rockchalk)
library(reghelper)
library(ggeffects)

apatheme <- theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line())

# ——————————————————————————————————————————————————————————————————————
# 0. Data preparation (same as Study1Analysis.R)
# ——————————————————————————————————————————————————————————————————————
S1.full <- read_csv("Study1Data.csv") %>%
  mutate(
    Date1 = as.Date(Date1, "%m/%d/%Y"),
    Date2 = as.Date(Date2, "%m/%d/%Y"),
    delay = Date2 - Date1,
    gender = factor(gender),
    amb = (neg + pos)/2 - abs(neg - pos),          # objective ambivalence
    attT1 = (att1 + att2 + att3)/3,                # T1 attitude
    attT2 = (att1_T2 + att2_T2 + att3_T2)/3,       # T2 attitude
    change = abs(attT1 - attT2),                   # attitude stability (absolute diff)
    samb = (subj1 + subj2 + subj3)/3,              # subjective ambivalence
    cert = select(., matches("cert[1-7]")) %>% rowMeans(na.rm = TRUE),
    ex = abs(attT1 - 4)                            # extremity
  )

S1 <- filter(S1.full, !is.na(attT2))

# ——————————————————————————————————————————————————————————————————————
# Part A: Standardized regression  (rjb_analysis.R, lines 115-129)
# ——————————————————————————————————————————————————————————————————————

cat("\n========== PART A: STANDARDIZED REGRESSION ==========\n")

S1 <- S1 %>%
  mutate(
    z_change    = c(scale(change)),
    z_cert      = c(scale(cert)),
    z_amb       = c(scale(amb)),
    z_ex        = c(scale(ex))
  )

# Step 1: main effects only
step1_z <- lm(z_change ~ z_cert + z_amb + z_ex, data = S1)
cat("\n--- Standardized Model — Main Effects ---\n")
print(summary(step1_z))

# Step 2: add interaction
step2_z <- lm(z_change ~ z_cert + z_amb + z_ex + z_cert:z_amb, data = S1)
cat("\n--- Standardized Model — With Interaction ---\n")
print(summary(step2_z))

# ANOVA comparison
cat("\n--- ANOVA: Step 1 vs Step 2 ---\n")
print(anova(step1_z, step2_z))

# Cohen's f^2 for the interaction
R2_1 <- summary(step1_z)$r.squared
R2_2 <- summary(step2_z)$r.squared
f2 <- (R2_2 - R2_1) / (1 - R2_2)
cat(sprintf("\nCohen's f^2 for interaction = %.4f\n", f2))


# --- Part A continued: Simple Slopes & Visualization ---
cat("\n\n--- Simple Slopes for Standardized Interaction (z_cert × z_amb) ---\n")

cat("\n--- rockchalk::plotSlopes ---\n")
plotSlopes(step2_z, modx = "z_amb", plotx = "z_cert", plotPoints = FALSE,
           n = 3, modxVals = "std.dev", col = c("black","black"), plotxRange = NULL)

cat("\n--- reghelper::simple_slopes ---\n")
print(simple_slopes(step2_z))

cat("\n--- ggplot2: simple slopes at ±1 SD amb ---\n")
plot_int_a <- ggpredict(step2_z, terms = c("z_cert [meansd]", "z_amb [meansd]")) %>%
  as.data.frame() %>%
  mutate(
    amb = factor(group, labels = c("Low Amb (-1 SD)", "Mean", "High Amb (+1 SD)")),
    cert = factor(x, labels = c("(-1 SD)", "Mean", "(+1 SD)"))
  ) %>%
  filter(amb != "Mean" & cert != "Mean") %>%
  mutate(amb = factor(amb))

pA <- ggplot(plot_int_a, aes(x = cert, y = predicted, group = amb, linetype = amb)) +
  geom_point(size = 3) +
  geom_line(size = 1) +
  scale_linetype_manual(values = c("dashed", "solid")) +
  labs(x = "Certainty (z)", y = "Attitude Change (z)", linetype = "") +
  apatheme + theme(legend.position = "top")
print(pA)


# ——————————————————————————————————————————————————————————————————————
# Part B: Three-way interaction  (LuttrellInferentialTestsOSF.R)
# ——————————————————————————————————————————————————————————————————————

cat("\n\n========== PART B: THREE-WAY INTERACTION ==========\n")

# Predict T2 attitude from T1 attitude * certainty * objective ambivalence
mod_3way_obj <- lm(attT2 ~ attT1 * cert * amb, data = S1)
cat("\n--- attT2 ~ attT1 * cert * amb (objective ambivalence) ---\n")
print(summary(mod_3way_obj))

# Predict T2 attitude from T1 attitude * certainty * subjective ambivalence
mod_3way_sub <- lm(attT2 ~ attT1 * cert * samb, data = S1)
cat("\n--- attT2 ~ attT1 * cert * samb (subjective ambivalence) ---\n")
print(summary(mod_3way_sub))

# ============================================================
# 三向交互可视化 (attT1 × cert × amb)
# 注意：ggpredict 返回的 x/group/facet 为数值型，不能按字符 "Low"/"High" 过滤
# ============================================================

cat("\n\n========== 三向交互可视化 (attT1 × cert × amb) ==========\n")

# --- 1. 客观矛盾性 (Objective Ambivalence) ---
plot_3way_obj <- ggpredict(mod_3way_obj,
                           terms = c("attT1 [meansd]", "cert [meansd]", "amb [meansd]")) %>%
  as.data.frame()

# 获取三个维度的唯一值（排序后: 1=-1SD, 2=Mean, 3=+1SD）
xv <- sort(unique(plot_3way_obj$x))
gv <- sort(unique(plot_3way_obj$group))
fv <- sort(unique(plot_3way_obj$facet))

plot_3way_obj <- plot_3way_obj %>%
  filter(x %in% xv[c(1, 3)], group %in% gv[c(1, 3)], facet %in% fv[c(1, 3)]) %>%
  mutate(
    attT1_level = factor(x,     levels = xv[c(1, 3)], labels = c("Low T1 (-1 SD)", "High T1 (+1 SD)")),
    cert_level  = factor(group, levels = gv[c(1, 3)], labels = c("Low Certainty (-1 SD)", "High Certainty (+1 SD)")),
    amb_level   = factor(facet, levels = fv[c(1, 3)], labels = c("Low Amb (-1 SD)", "High Amb (+1 SD)"))
  )

p_3way_obj <- ggplot(plot_3way_obj,
                     aes(x = attT1_level, y = predicted,
                         group = cert_level, linetype = cert_level)) +
  geom_point(size = 3) +
  geom_line(size = 1) +
  facet_wrap(~ amb_level) +
  scale_linetype_manual(values = c("dashed", "solid")) +
  labs(x = "T1 Attitude", y = "Predicted T2 Attitude",
       linetype = "Certainty", title = "Objective Ambivalence") +
  apatheme +
  theme(legend.position = "top",
        strip.background = element_blank(),
        strip.text = element_text(size = 10, face = "bold"))
print(p_3way_obj)


# --- 2. 主观矛盾性 (Subjective Ambivalence) ---
plot_3way_sub <- ggpredict(mod_3way_sub,
                           terms = c("attT1 [meansd]", "cert [meansd]", "samb [meansd]")) %>%
  as.data.frame()

xv2 <- sort(unique(plot_3way_sub$x))
gv2 <- sort(unique(plot_3way_sub$group))
fv2 <- sort(unique(plot_3way_sub$facet))

plot_3way_sub <- plot_3way_sub %>%
  filter(x %in% xv2[c(1, 3)], group %in% gv2[c(1, 3)], facet %in% fv2[c(1, 3)]) %>%
  mutate(
    attT1_level = factor(x,     levels = xv2[c(1, 3)], labels = c("Low T1 (-1 SD)", "High T1 (+1 SD)")),
    cert_level  = factor(group, levels = gv2[c(1, 3)], labels = c("Low Certainty (-1 SD)", "High Certainty (+1 SD)")),
    samb_level  = factor(facet, levels = fv2[c(1, 3)], labels = c("Low Subj. Amb (-1 SD)", "High Subj. Amb (+1 SD)"))
  )

p_3way_sub <- ggplot(plot_3way_sub,
                     aes(x = attT1_level, y = predicted,
                         group = cert_level, linetype = cert_level)) +
  geom_point(size = 3) +
  geom_line(size = 1) +
  facet_wrap(~ samb_level) +
  scale_linetype_manual(values = c("dashed", "solid")) +
  labs(x = "T1 Attitude", y = "Predicted T2 Attitude",
       linetype = "Certainty", title = "Subjective Ambivalence") +
  apatheme +
  theme(legend.position = "top",
        strip.background = element_blank(),
        strip.text = element_text(size = 10, face = "bold"))
print(p_3way_sub)


# ——————————————————————————————————————————————————————————————————————
# Part C: Bayesian estimation
# ——————————————————————————————————————————————————————————————————————

cat("\n\n========== PART C: BAYESIAN ESTIMATION ==========\n")
cat(" Bayes Factor (BF10): comparing interaction model vs. main-effects-only model\n")
cat(" Posterior distribution & 95% HDI for the interaction coefficient\n")
cat(" Sensitivity analysis across prior scales\n\n")

# ——————————————————————————————————————————————————————————————————————
# C1. Bayes Factor model comparison
# ——————————————————————————————————————————————————————————————————————

# Remove rows with missing values for the variables used
bf_data <- S1 %>% drop_na(change, ex, amb, cert)

# Main-effects model
bf_main <- lmBF(change ~ ex + amb + cert, data = as.data.frame(bf_data))
# Interaction model
bf_int  <- lmBF(change ~ ex + amb + cert + amb:cert, data = as.data.frame(bf_data))

# BF10: evidence for interaction model over main-effects model
BF10 <- bf_int / bf_main
bf10_val <- as.vector(BF10)
cat(sprintf("\nBayes Factor (BF10, interaction vs. main effects) = %.4f\n", bf10_val))
cat(sprintf("Interpretation: The data are %.2f times more likely under the interaction model\n", bf10_val))
cat(sprintf("               than under the main-effects-only model.\n"))
if (bf10_val > 3) {
  cat("               -> Substantial evidence FOR the interaction.\n")
} else if (bf10_val > 1) {
  cat("               -> Anecdotal evidence FOR the interaction.\n")
} else if (bf10_val > 1/3) {
  cat("               -> Anecdotal evidence AGAINST the interaction.\n")
} else {
  cat("               -> Substantial evidence AGAINST the interaction.\n")
}

# Also compute BF01 (the reciprocal)
cat(sprintf("BF01 (main effects / interaction) = %.4f\n", 1 / bf10_val))


# ——————————————————————————————————————————————————————————————————————
# C2. Posterior distribution & 95% HDI for the interaction coefficient
# ——————————————————————————————————————————————————————————————————————

set.seed(42)
chains <- posterior(bf_int, iterations = 10000, progress = FALSE)

# Column names in the BayesFactor posterior
cn <- colnames(chains)
cn_clean <- cn[!is.na(cn)]
cat("Posterior column names:", paste(cn_clean, collapse = ", "), "\n")

int_col <- grep("amb.*cert|cert.*amb", cn_clean, value = TRUE)
if (length(int_col) == 0) {
  int_col <- grep("amb", cn_clean, value = TRUE)
  int_col <- grep("cert", int_col, value = TRUE)
}
cat(sprintf("\nInteraction coefficient column: %s\n", int_col))

int_samples <- chains[, int_col]

# Summary statistics
post_mean <- mean(int_samples)
post_med  <- median(int_samples)
post_sd   <- sd(int_samples)

# 95% HDI
hdi_result <- hdi(int_samples, ci = 0.95)
hdi_lower <- hdi_result$CI_low
hdi_upper <- hdi_result$CI_high

# Does zero fall outside the HDI?
zero_in_hdi <- (hdi_lower < 0 & hdi_upper > 0)

cat(sprintf("\nPosterior of interaction coefficient (amb:cert):\n"))
cat(sprintf("  Mean     = %.4f\n", post_mean))
cat(sprintf("  Median   = %.4f\n", post_med))
cat(sprintf("  SD       = %.4f\n", post_sd))
cat(sprintf("  95%% HDI  = [%.4f, %.4f]\n", hdi_lower, hdi_upper))
cat(sprintf("  Zero %s inside the 95%% HDI\n",
            ifelse(zero_in_hdi, "IS", "is NOT")))

# Probability that the interaction coefficient is negative/positive
p_neg <- mean(int_samples < 0)
p_pos <- mean(int_samples > 0)
cat(sprintf("  P(coefficient < 0) = %.3f\n", p_neg))
cat(sprintf("  P(coefficient > 0) = %.3f\n", p_pos))


# ——————————————————————————————————————————————————————————————————————
# C3. Sensitivity analysis: vary prior scale (rscale)
# ——————————————————————————————————————————————————————————————————————

cat("\n\n--- Sensitivity Analysis ---\n")
cat("Varying the JZS prior scale parameter (rscale):\n")
cat("  wide   = 1.00 (diffuse)\n")
cat("  medium = 0.50 (default)\n")
cat("  narrow = 0.20 (conservative)\n\n")

rscale_values <- c(0.20, 0.50, 1.00)
rscale_names  <- c("narrow (0.20)", "medium (0.50)", "wide (1.00)")

sensitivity_table <- data.frame(
  Prior_Scale = character(),
  BF10        = numeric(),
  Post_Mean   = numeric(),
  HDI_Lower   = numeric(),
  HDI_Upper   = numeric(),
  stringsAsFactors = FALSE
)

for (i in seq_along(rscale_values)) {
  r <- rscale_values[i]

  bf_main_r <- lmBF(change ~ ex + amb + cert,
                     data = as.data.frame(bf_data),
                     rscaleFixed = r)
  bf_int_r  <- lmBF(change ~ ex + amb + cert + amb:cert,
                     data = as.data.frame(bf_data),
                     rscaleFixed = r)

  bf10_r <- as.vector(bf_int_r / bf_main_r)

  # Posterior under this prior
  set.seed(42)
  chains_r <- posterior(bf_int_r, iterations = 10000, progress = FALSE)
  cn_r <- colnames(chains_r)
  cn_r <- cn_r[!is.na(cn_r)]
  int_col_r <- grep("amb.*cert|cert.*amb", cn_r, value = TRUE)
  if (length(int_col_r) == 0) {
    int_col_r <- grep("amb", cn_r, value = TRUE)
    int_col_r <- grep("cert", int_col_r, value = TRUE)
  }
  int_samples_r <- chains_r[, int_col_r]
  post_mean_r <- mean(int_samples_r)
  hdi_r <- hdi(int_samples_r, ci = 0.95)

  sensitivity_table <- rbind(sensitivity_table, data.frame(
    Prior_Scale = rscale_names[i],
    BF10        = round(bf10_r, 4),
    Post_Mean   = round(post_mean_r, 4),
    HDI_Lower   = round(hdi_r$CI_low, 4),
    HDI_Upper   = round(hdi_r$CI_high, 4),
    stringsAsFactors = FALSE
  ))
}

print(sensitivity_table, row.names = FALSE)

cat("\nConclusion across priors:\n")
if (all(sensitivity_table$BF10 > 1)) {
  cat("  BF10 > 1 across all prior scales — the data consistently favor\n")
  cat("  the interaction model over the main-effects-only model.\n")
} else {
  cat("  BF10 crosses 1 across prior scales — the evidence is sensitive\n")
  cat("  to prior specification; caution is warranted.\n")
}
# Check if HDI results are available (not NA)
hdi_available <- !any(is.na(sensitivity_table$HDI_Lower))
if (hdi_available) {
  if (all(sensitivity_table$HDI_Lower < 0 & sensitivity_table$HDI_Upper > 0)) {
    cat("  The 95% HDI contains zero across all priors — the interaction\n")
    cat("  effect is not precisely estimated.\n")
  } else {
    cat("  The 95% HDI excludes zero under at least one prior — suggesting\n")
    cat("  a credible non-zero interaction effect.\n")
  }
}


# ——————————————————————————————————————————————————————————————————————
# Part D: Linear Mixed Model (LMM) — wide-to-long, Time × Cert × Amb
# ——————————————————————————————————————————————————————————————————————

cat("\n\n========== PART D: LINEAR MIXED MODEL ==========\n")
cat(" Reshape data to long format; use raw attitude scores (not absolute difference)\n")
cat(" Time × Certainty × Ambivalence three-way interaction\n")
cat(" Random intercept + random slope for Time per participant\n")
cat(" Using ALL available data (including T1-only)\n\n")

library(lme4)
library(lmerTest)

# ——————————————————————————————————————————————————————————————
# D0. Prepare long-format data
# ——————————————————————————————————————————————————————————————

# Center continuous predictors at grand mean for interpretability
S1.full <- S1.full %>%
  mutate(
    cert_c = cert - mean(cert, na.rm = TRUE),
    amb_c  = amb  - mean(amb,  na.rm = TRUE),
    ex_c   = ex   - mean(ex,   na.rm = TRUE)
  )

# Reshape: each row = one time-point per participant
long_data <- S1.full %>%
  select(ID, cert_c, amb_c, ex_c, att1, att2, att3, att1_T2, att2_T2, att3_T2) %>%
  pivot_longer(
    cols      = c(att1, att2, att3, att1_T2, att2_T2, att3_T2),
    names_to  = "item",
    values_to = "attitude"
  ) %>%
  mutate(
    time = ifelse(str_detect(item, "_T2$"), 1, 0),   # 0 = T1, 1 = T2
    # contrast code: -0.5 = T1, +0.5 = T2 (so intercept = grand mean, slope = change)
    time_c = time - 0.5
  ) %>%
  filter(!is.na(attitude))   # drop missing observations

cat(sprintf("Long-format data: %d observations from %d participants\n",
            nrow(long_data), n_distinct(long_data$ID)))
cat(sprintf("  T1 observations: %d\n", sum(long_data$time == 0)))
cat(sprintf("  T2 observations: %d\n", sum(long_data$time == 1)))
cat(sprintf("  Participants with both T1 & T2: %d\n",
            long_data %>% group_by(ID) %>% summarise(n = n_distinct(time)) %>% filter(n == 2) %>% nrow()))
cat(sprintf("  Participants with T1 only: %d\n",
            long_data %>% group_by(ID) %>% summarise(n = n_distinct(time)) %>% filter(n == 1) %>% nrow()))


# ——————————————————————————————————————————————————————————————
# D1. Model building
# ——————————————————————————————————————————————————————————————

# M0: Null (intercept only)
cat("\n--- Model D0: Random intercept only ---\n")
m0 <- lmer(attitude ~ 1 + (1 | ID), data = long_data, REML = FALSE)
print(summary(m0, correlation = FALSE))

# M1: Add fixed effects of time, cert, amb, ex
cat("\n--- Model D1: Main effects (time + cert + amb + ex) ---\n")
m1 <- lmer(attitude ~ time_c + cert_c + amb_c + ex_c + (1 | ID),
           data = long_data, REML = FALSE)
print(summary(m1, correlation = FALSE))

# M2: Add two-way interactions
cat("\n--- Model D2: Two-way interactions ---\n")
m2 <- lmer(attitude ~ time_c * cert_c + time_c * amb_c + cert_c * amb_c + ex_c +
           (1 | ID),
           data = long_data, REML = FALSE)
print(summary(m2, correlation = FALSE))

# M3: Full three-way interaction
cat("\n--- Model D3: Time × Cert × Amb (three-way) ---\n")
m3 <- lmer(attitude ~ time_c * cert_c * amb_c + ex_c + (1 | ID),
           data = long_data, REML = FALSE)
print(summary(m3, correlation = FALSE))


# ——————————————————————————————————————————————————————————————
# D2. Model comparison (likelihood ratio tests)
# ——————————————————————————————————————————————————————————————

cat("\n--- Likelihood Ratio Tests ---\n")

cat("\nD1 vs D0 (main effects vs null):\n")
print(anova(m0, m1))

cat("\nD2 vs D1 (adding two-way interactions):\n")
print(anova(m1, m2))

cat("\nD3 vs D2 (adding three-way interaction):\n")
print(anova(m2, m3))


# ——————————————————————————————————————————————————————————————
# D3. Add random slope for time
# ——————————————————————————————————————————————————————————————

cat("\n--- Model D4: With random slope for time (uncorrelated) ---\n")
m4 <- lmer(attitude ~ time_c * cert_c * amb_c + ex_c +
           (1 | ID) + (0 + time_c | ID),
           data = long_data, REML = FALSE)
print(summary(m4, correlation = FALSE))

cat("\nD4 vs D3 (random slope test):\n")
print(anova(m3, m4))

cat("\n--- Model D5: With correlated random slope for time ---\n")
m5 <- lmer(attitude ~ time_c * cert_c * amb_c + ex_c +
           (1 + time_c | ID),
           data = long_data, REML = FALSE)
print(summary(m5, correlation = FALSE))

cat("\nD5 vs D3 (correlated random slope test):\n")
print(anova(m3, m5))


# ——————————————————————————————————————————————————————————————
# D4. Best model summary (with REML for unbiased estimates)
# ——————————————————————————————————————————————————————————————

# Refit best model with REML = TRUE
cat("\n--- Final Model (REML): attitude ~ time_c * cert_c * amb_c + ex_c + (1 | ID) ---\n")
m_final <- lmer(attitude ~ time_c * cert_c * amb_c + ex_c + (1 | ID),
                data = long_data, REML = TRUE)
print(summary(m_final, correlation = FALSE))

# Extract the three-way interaction
threeway_fixef <- summary(m_final)$coefficients
threeway_row <- grep("time_c:cert_c:amb_c|cert_c:time_c:amb_c", rownames(threeway_fixef))
cat("\n--- Three-way interaction (Time × Cert × Amb) ---\n")
if (length(threeway_row) > 0) {
  est <- threeway_fixef[threeway_row, "Estimate"]
  se  <- threeway_fixef[threeway_row, "Std. Error"]
  t   <- threeway_fixef[threeway_row, "t value"]
  p   <- threeway_fixef[threeway_row, "Pr(>|t|)"]
  cat(sprintf("  Estimate = %.4f, SE = %.4f, t = %.3f, p = %.4f\n", est, se, t, p))
  cat(sprintf("  95%% CI = [%.4f, %.4f]\n", est - 1.96*se, est + 1.96*se))
}

# Also report the critical Time × Amb interaction (if the three-way is the focal)
# and the Time × Cert interaction
time_amb_row <- grep("time_c:amb_c", rownames(threeway_fixef))
time_amb_row <- time_amb_row[!grepl("cert", names(time_amb_row))]
cat("\n--- Time × Amb interaction (averaged over cert) ---\n")
if (length(time_amb_row) > 0) {
  est <- threeway_fixef[time_amb_row[1], "Estimate"]
  se  <- threeway_fixef[time_amb_row[1], "Std. Error"]
  t   <- threeway_fixef[time_amb_row[1], "t value"]
  p   <- threeway_fixef[time_amb_row[1], "Pr(>|t|)"]
  cat(sprintf("  Estimate = %.4f, SE = %.4f, t = %.3f, p = %.4f\n", est, se, t, p))
}

# ——————————————————————————————————————————————————————————————
# D5. Variance components (random effects)
# ——————————————————————————————————————————————————————————————

cat("\n--- Variance Components ---\n")
vc <- VarCorr(m_final)
print(vc, comp = c("Variance", "Std.Dev."))

# ICC
icc_val <- as.numeric(vc$ID) / (as.numeric(vc$ID) + attr(vc, "sc")^2)
cat(sprintf("ICC (proportion of variance at person level) = %.4f\n", icc_val))


# ——————————————————————————————————————————————————————————————
# D6. Simple slopes / conditional effects decomposition
# ——————————————————————————————————————————————————————————————

cat("\n--- Conditional effect of time at ±1 SD of cert and amb ---\n")

cert_sd <- sd(S1.full$cert, na.rm = TRUE)
amb_sd  <- sd(S1.full$amb,  na.rm = TRUE)
cert_m  <- mean(S1.full$cert, na.rm = TRUE)
amb_m   <- mean(S1.full$amb,  na.rm = TRUE)

# Create grid for predictions
cond_grid <- expand.grid(
  time_c = c(-0.5, 0.5),
  cert_c = c(-cert_sd, 0, cert_sd),
  amb_c  = c(-amb_sd, 0, amb_sd),
  ex_c   = 0
)

# Get predictions from the final model
cond_grid$pred <- predict(m_final, newdata = cond_grid, re.form = NA)

# Compute conditional Time slopes (T2 - T1) for each cert × amb combination
cond_slopes <- cond_grid %>%
  group_by(cert_c, amb_c) %>%
  summarise(
    time_slope = pred[time_c == 0.5] - pred[time_c == -0.5],
    .groups = "drop"
  ) %>%
  mutate(
    cert_level = case_when(
      cert_c == -cert_sd ~ "-1 SD",
      cert_c == 0        ~ "Mean",
      cert_c == cert_sd  ~ "+1 SD"
    ),
    amb_level = case_when(
      amb_c == -amb_sd  ~ "-1 SD",
      amb_c == 0        ~ "Mean",
      amb_c == amb_sd   ~ "+1 SD"
    )
  )

cat("\nConditional Time slopes (attitude change from T1 to T2):\n")
print(as.data.frame(cond_slopes[, c("cert_level", "amb_level", "time_slope")]))


cat("\n========== ALL ANALYSES COMPLETE ==========\n")
