########################################################################
# Luttrell, Petty, & Brinol (2016) Study 2 — Other Analyses
# Includes:
#   Part A: Standardized regression
#   Part B: Three-way interaction
#   Part C: Bayesian estimation (BF10, HDI, sensitivity)
#   Part D: Linear Mixed Model (LMM) with 3 time points
########################################################################

library(tidyverse)
library(psych)
library(Hmisc)
library(BayesFactor)
library(bayestestR)
library(lme4)
library(lmerTest)

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
# 0. Data preparation (same as Study2Analysis.R)
# ——————————————————————————————————————————————————————————————————————
S2.full <- read_csv("Study2Data.csv") %>%
  mutate(
    Date1 = as.Date(Date1, "%m/%d/%Y"),
    Date2 = as.Date(Date2, "%m/%d/%Y"),
    Date3 = as.Date(Date3, "%m/%d/%Y"),
    delay.12 = Date2 - Date1,
    delay = Date3 - Date1,
    gender = factor(gender),
    amb = (neg + pos)/2 - abs(neg - pos),
    attT1 = (att1 + att2 + att3)/3,
    attT2 = (att1_T2 + att2_T2 + att3_T2)/3,
    attT3 = (att1_T3 + att2_T3 + att3_T3)/3,
    change = abs(attT3 - attT1),         # paper focus: one-year delay
    change.12 = abs(attT1 - attT2),       # shorter 2-month delay
    samb = (subj1 + subj2 + subj3)/3,
    cert = select(., matches("cert[1-7]")) %>% rowMeans(na.rm = TRUE),
    ex = abs(attT1 - 4)
  )

S2 <- filter(S2.full, !is.na(attT3))

# ——————————————————————————————————————————————————————————————————————
# Part A: Standardized regression
# ——————————————————————————————————————————————————————————————————————

cat("\n========== PART A: STANDARDIZED REGRESSION ==========\n")

S2 <- S2 %>%
  mutate(
    z_change    = c(scale(change)),
    z_cert      = c(scale(cert)),
    z_amb       = c(scale(amb)),
    z_ex        = c(scale(ex))
  )

step1_z <- lm(z_change ~ z_cert + z_amb + z_ex, data = S2)
cat("\n--- Step 1: Main effects (one-year change) ---\n")
print(summary(step1_z))

step2_z <- lm(z_change ~ z_cert + z_amb + z_ex + z_cert:z_amb, data = S2)
cat("\n--- Step 2: + Interaction ---\n")
print(summary(step2_z))

cat("\n--- ANOVA ---\n")
print(anova(step1_z, step2_z))

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


# Also do 2-month change for comparison
cat("\n\n--- Two-month change (change.12) ---\n")
S2.12 <- S2 %>% filter(!is.na(change.12)) %>%
  mutate(
    z_change12 = c(scale(change.12)),
    z_cert12   = c(scale(cert)),
    z_amb12    = c(scale(amb)),
    z_ex12     = c(scale(ex))
  )
step2_z12 <- lm(z_change12 ~ z_cert12 + z_amb12 + z_ex12 + z_cert12:z_amb12,
                data = S2.12)
print(summary(step2_z12))

# ——————————————————————————————————————————————————————————————————————
# Part B: Three-way interaction
# ——————————————————————————————————————————————————————————————————————

cat("\n\n========== PART B: THREE-WAY INTERACTION ==========\n")

# Predict T3 attitude from T1 attitude * certainty * ambivalence
mod_3way_obj <- lm(attT3 ~ attT1 * cert * amb, data = S2)
cat("\n--- attT3 ~ attT1 * cert * amb (objective amb) ---\n")
print(summary(mod_3way_obj))

mod_3way_sub <- lm(attT3 ~ attT1 * cert * samb, data = S2)
cat("\n--- attT3 ~ attT1 * cert * samb (subjective amb) ---\n")
print(summary(mod_3way_sub))


# --- Part B.4: Three-way Interaction Visualization (attT1 × cert × amb) ---
cat("\n\n========== Part B.4: THREE-WAY PLOTS — attT1 × cert × amb ==========\n")

# --- 1. 客观矛盾性 (Objective Ambivalence) ---
plot_3way_obj <- ggpredict(mod_3way_obj,
                           terms = c("attT1 [meansd]", "cert [meansd]", "amb [meansd]")) %>%
  as.data.frame()

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
  labs(x = "T1 Attitude", y = "Predicted T3 Attitude",
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
  labs(x = "T1 Attitude", y = "Predicted T3 Attitude",
       linetype = "Certainty", title = "Subjective Ambivalence") +
  apatheme +
  theme(legend.position = "top",
        strip.background = element_blank(),
        strip.text = element_text(size = 10, face = "bold"))
print(p_3way_sub)


# --- Part B.5: VIF (Multicollinearity) & Bootstrap for Three-way Interaction ---
cat("\n\n========== Part B.5: VIF & BOOTSTRAP ==========\n")
library(car)

cat("\n--- VIF for attT3 ~ attT1 * cert * amb (objective, three-way) ---\n")
print(vif(mod_3way_obj))

cat("\n--- VIF for attT3 ~ attT1 * cert * samb (subjective, three-way) ---\n")
print(vif(mod_3way_sub))

cat("\nNote: High VIF for interaction terms is expected due to product-term\n")
cat("collinearity with their constituents; mean-centering predictors helps.\n")

# --- Bootstrap 95% CI (5000 resamples) for three-way interaction ---
cat("\n\n--- Bootstrap (5000 reps) — Three-way interaction coefficients ---\n")
library(boot)

# Bootstrap function: objective ambivalence three-way
boot_3way_obj <- function(data, indices) {
  d <- data[indices, ]
  m <- lm(attT3 ~ attT1 * cert * amb, data = d)
  coef(m)["attT1:cert:amb"]
}

set.seed(42)
boot_obj <- boot(S2, boot_3way_obj, R = 5000)
ci_obj <- boot.ci(boot_obj, type = c("norm", "perc", "bca"))

cat(sprintf("\nObjective Amb — three-way attT1:cert:amb:\n"))
cat(sprintf("  Original coef = %.4f\n", boot_obj$t0))
cat(sprintf("  Bootstrap SE  = %.4f\n", sd(boot_obj$t)))
cat(sprintf("  95%% CI (normal):      [%.4f, %.4f]\n",
            ci_obj$normal[2], ci_obj$normal[3]))
cat(sprintf("  95%% CI (percentile):  [%.4f, %.4f]\n",
            ci_obj$percent[4], ci_obj$percent[5]))
cat(sprintf("  95%% CI (BCa):         [%.4f, %.4f]\n",
            ci_obj$bca[4], ci_obj$bca[5]))

# Bootstrap function: subjective ambivalence three-way
boot_3way_sub <- function(data, indices) {
  d <- data[indices, ]
  m <- lm(attT3 ~ attT1 * cert * samb, data = d)
  coef(m)["attT1:cert:samb"]
}

set.seed(42)
boot_sub <- boot(S2, boot_3way_sub, R = 5000)
ci_sub <- boot.ci(boot_sub, type = c("norm", "perc", "bca"))

cat(sprintf("\nSubjective Amb — three-way attT1:cert:samb:\n"))
cat(sprintf("  Original coef = %.4f\n", boot_sub$t0))
cat(sprintf("  Bootstrap SE  = %.4f\n", sd(boot_sub$t)))
cat(sprintf("  95%% CI (normal):      [%.4f, %.4f]\n",
            ci_sub$normal[2], ci_sub$normal[3]))
cat(sprintf("  95%% CI (percentile):  [%.4f, %.4f]\n",
            ci_sub$percent[4], ci_sub$percent[5]))
cat(sprintf("  95%% CI (BCa):         [%.4f, %.4f]\n",
            ci_sub$bca[4], ci_sub$bca[5]))



# ——————————————————————————————————————————————————————————————————————
# Part C: Bayesian estimation
# ——————————————————————————————————————————————————————————————————————

cat("\n\n========== PART C: BAYESIAN ESTIMATION ==========\n")

# ——————————————————————————————————————————————————————————————
# C1. Bayes Factor
# ——————————————————————————————————————————————————————————————

bf_data <- S2 %>% drop_na(change, ex, amb, cert)

bf_main <- lmBF(change ~ ex + amb + cert, data = as.data.frame(bf_data))
bf_int  <- lmBF(change ~ ex + amb + cert + amb:cert, data = as.data.frame(bf_data))

BF10 <- bf_int / bf_main
bf10_val <- as.vector(BF10)
cat(sprintf("\nBF10 (interaction vs. main effects) = %.4f\n", bf10_val))
cat(sprintf("BF01 = %.4f\n", 1 / bf10_val))
if (bf10_val > 3) {
  cat("-> Substantial evidence FOR the interaction.\n")
} else if (bf10_val > 1) {
  cat("-> Anecdotal evidence FOR the interaction.\n")
} else if (bf10_val > 1/3) {
  cat("-> Anecdotal evidence AGAINST the interaction.\n")
} else {
  cat("-> Substantial evidence AGAINST the interaction.\n")
}

# ——————————————————————————————————————————————————————————————
# C2. Posterior & HDI
# ——————————————————————————————————————————————————————————————

set.seed(42)
chains <- posterior(bf_int, iterations = 10000, progress = FALSE)

cn_clean <- colnames(chains); cn_clean <- cn_clean[!is.na(cn_clean)]
cat("Posterior columns:", paste(cn_clean, collapse = ", "), "\n")

int_col <- grep("amb.*cert|cert.*amb", cn_clean, value = TRUE)
if (length(int_col) == 0) {
  int_col <- grep("amb", cn_clean, value = TRUE)
  int_col <- grep("cert", int_col, value = TRUE)
}
int_samples <- chains[, int_col]

post_mean <- mean(int_samples)
post_med  <- median(int_samples)
post_sd   <- sd(int_samples)
hdi_res   <- hdi(int_samples, ci = 0.95)
hdi_lower <- hdi_res$CI_low
hdi_upper <- hdi_res$CI_high
zero_in   <- hdi_lower < 0 & hdi_upper > 0

cat(sprintf("\nInteraction coefficient (amb:cert):\n"))
cat(sprintf("  Mean = %.4f, Median = %.4f, SD = %.4f\n", post_mean, post_med, post_sd))
cat(sprintf("  95%% HDI = [%.4f, %.4f]\n", hdi_lower, hdi_upper))
cat(sprintf("  Zero %s inside HDI\n", ifelse(zero_in, "IS", "is NOT")))
cat(sprintf("  P(coef < 0) = %.3f, P(coef > 0) = %.3f\n",
            mean(int_samples < 0), mean(int_samples > 0)))

# ——————————————————————————————————————————————————————————————
# C3. Sensitivity analysis
# ——————————————————————————————————————————————————————————————

cat("\n--- Sensitivity Analysis ---\n")
rscale_vals <- c(0.20, 0.50, 1.00)
rscale_names <- c("narrow (0.20)", "medium (0.50)", "wide (1.00)")

sens <- data.frame(Prior_Scale = character(), BF10 = numeric(),
                   Post_Mean = numeric(), HDI_Lower = numeric(),
                   HDI_Upper = numeric(), stringsAsFactors = FALSE)

for (i in seq_along(rscale_vals)) {
  r <- rscale_vals[i]
  bf_m <- lmBF(change ~ ex + amb + cert, data = as.data.frame(bf_data), rscaleFixed = r)
  bf_i <- lmBF(change ~ ex + amb + cert + amb:cert, data = as.data.frame(bf_data), rscaleFixed = r)
  bf10_r <- as.vector(bf_i / bf_m)

  set.seed(42)
  ch_r <- posterior(bf_i, iterations = 10000, progress = FALSE)
  cn_r <- colnames(ch_r); cn_r <- cn_r[!is.na(cn_r)]
  ic_r <- grep("amb.*cert|cert.*amb", cn_r, value = TRUE)
  if (length(ic_r) == 0) {
    ic_r <- grep("amb", cn_r, value = TRUE); ic_r <- grep("cert", ic_r, value = TRUE)
  }
  s_r <- ch_r[, ic_r]
  hdi_r <- hdi(s_r, ci = 0.95)

  sens <- rbind(sens, data.frame(
    Prior_Scale = rscale_names[i], BF10 = round(bf10_r, 4),
    Post_Mean = round(mean(s_r), 4),
    HDI_Lower = round(hdi_r$CI_low, 4), HDI_Upper = round(hdi_r$CI_high, 4),
    stringsAsFactors = FALSE
  ))
}

print(sens, row.names = FALSE)

if (all(sens$BF10 > 1)) {
  cat("BF10 > 1 across all priors — consistent evidence for interaction.\n")
} else {
  cat("BF10 crosses 1 — evidence is sensitive to prior choice.\n")
}

# ——————————————————————————————————————————————————————————————————————
# Part D: Linear Mixed Model (3 time points!)
# ——————————————————————————————————————————————————————————————————————

cat("\n\n========== PART D: LINEAR MIXED MODEL ==========\n")
cat(" Long-format with 3 time points (T1, T2, T3)\n")
cat(" Random intercept + random slope per participant\n")
cat(" Time × Certainty × Ambivalence interaction\n")
cat(" Using ALL available data\n\n")

# ——————————————————————————————————————————————————————————————
# D0. Prepare long data
# ——————————————————————————————————————————————————————————————

S2.full <- S2.full %>%
  mutate(
    cert_c = cert - mean(cert, na.rm = TRUE),
    amb_c  = amb  - mean(amb,  na.rm = TRUE),
    ex_c   = ex   - mean(ex,   na.rm = TRUE)
  )

long_data <- S2.full %>%
  select(ID, cert_c, amb_c, ex_c,
         att1, att2, att3,
         att1_T2, att2_T2, att3_T2,
         att1_T3, att2_T3, att3_T3) %>%
  pivot_longer(
    cols      = c(att1, att2, att3, att1_T2, att2_T2, att3_T2, att1_T3, att2_T3, att3_T3),
    names_to  = "item",
    values_to = "attitude"
  ) %>%
  mutate(
    time = case_when(
      str_detect(item, "_T3$") ~ 2,
      str_detect(item, "_T2$") ~ 1,
      TRUE                     ~ 0
    ),
    # 实际经过月份: T1=0月(基线), T2=2月, T3=12月(一年延迟)
    # 中心化使截距代表跨时间点的态度均值
    time_month = case_when(
      time == 0 ~ 0,
      time == 1 ~ 2,
      time == 2 ~ 12
    )
  ) %>%
  filter(!is.na(attitude)) %>%
  mutate(time_c = time_month - mean(time_month))

cat(sprintf("Long data: %d obs from %d participants\n",
            nrow(long_data), n_distinct(long_data$ID)))
cat(sprintf("  T1: %d, T2: %d, T3: %d\n",
            sum(long_data$time == 0), sum(long_data$time == 1), sum(long_data$time == 2)))
cat(sprintf("  Time_c range: [%.2f, %.2f] (months, centered)\n",
            min(long_data$time_c), max(long_data$time_c)))
cat(sprintf("  Participants with T1+T2+T3: %d\n",
            long_data %>% group_by(ID) %>% summarise(n = n_distinct(time)) %>% filter(n == 3) %>% nrow()))
cat(sprintf("  Participants with T1 only: %d\n",
            long_data %>% group_by(ID) %>% summarise(n = n_distinct(time)) %>% filter(n == 1) %>% nrow()))

# ——————————————————————————————————————————————————————————————
# D1. Model building
# ——————————————————————————————————————————————————————————————

m0 <- lmer(attitude ~ 1 + (1 | ID), data = long_data, REML = FALSE)
cat("\n--- M0: Null (random intercept) ---\n")
cat(sprintf("AIC = %.1f, logLik = %.1f\n", AIC(m0), logLik(m0)))

m1 <- lmer(attitude ~ time_c + cert_c + amb_c + ex_c + (1 | ID),
           data = long_data, REML = FALSE)
cat("\n--- M1: Main effects ---\n")
print(summary(m1, correlation = FALSE))

m2 <- lmer(attitude ~ time_c * cert_c + time_c * amb_c + cert_c * amb_c + ex_c + (1 | ID),
           data = long_data, REML = FALSE)
cat("\n--- M2: Two-way interactions ---\n")
print(summary(m2, correlation = FALSE))

m3 <- lmer(attitude ~ time_c * cert_c * amb_c + ex_c + (1 | ID),
           data = long_data, REML = FALSE)
cat("\n--- M3: Three-way interaction ---\n")
print(summary(m3, correlation = FALSE))

# ——————————————————————————————————————————————————————————————
# D2. Likelihood ratio tests
# ——————————————————————————————————————————————————————————————

cat("\n--- Likelihood Ratio Tests ---\n")
cat("\nM1 vs M0 (main effects vs null):\n")
print(anova(m0, m1))
cat("\nM2 vs M1 (two-way interactions):\n")
print(anova(m1, m2))
cat("\nM3 vs M2 (three-way):\n")
print(anova(m2, m3))

# ——————————————————————————————————————————————————————————————
# D3. Random slope for time
# ——————————————————————————————————————————————————————————————

cat("\n--- M4: With uncorrelated random slope for time ---\n")
m4 <- lmer(attitude ~ time_c * cert_c * amb_c + ex_c +
           (1 | ID) + (0 + time_c | ID),
           data = long_data, REML = FALSE)
cat(sprintf("AIC = %.1f, logLik = %.1f\n", AIC(m4), logLik(m4)))
cat("Random effects:\n")
print(VarCorr(m4), comp = c("Variance", "Std.Dev."))
cat("\nFixed effects:\n")
print(coef(summary(m4)))

cat("\nM4 vs M3 (random slope test):\n")
print(anova(m3, m4))

cat("\n--- M5: Correlated random slope ---\n")
m5 <- lmer(attitude ~ time_c * cert_c * amb_c + ex_c +
           (1 + time_c | ID),
           data = long_data, REML = FALSE)
cat(sprintf("AIC = %.1f, logLik = %.1f\n", AIC(m5), logLik(m5)))
print(summary(m5, correlation = FALSE))

# ——————————————————————————————————————————————————————————————
# D4. Final REML model
# ——————————————————————————————————————————————————————————————

cat("\n--- Final Model (REML): (1 | ID) ---\n")
m_final <- lmer(attitude ~ time_c * cert_c * amb_c + ex_c + (1 | ID),
                data = long_data, REML = TRUE)
print(summary(m_final, correlation = FALSE))

# Three-way interaction
fixef_tab <- summary(m_final)$coefficients
tw_row <- grep("time_c:cert_c:amb_c|cert_c:time_c:amb_c", rownames(fixef_tab))
cat("\n--- Three-way interaction (Time × Cert × Amb) ---\n")
if (length(tw_row) > 0) {
  e <- fixef_tab[tw_row, "Estimate"]
  s <- fixef_tab[tw_row, "Std. Error"]
  t <- fixef_tab[tw_row, "t value"]
  p <- fixef_tab[tw_row, "Pr(>|t|)"]
  cat(sprintf("  Est = %.4f, SE = %.4f, t = %.3f, p = %.4f\n", e, s, t, p))
  cat(sprintf("  95%% CI = [%.4f, %.4f]\n", e - 1.96*s, e + 1.96*s))
}

# Time × Amb interaction
ta_row <- grep("^time_c:amb_c$", rownames(fixef_tab))
cat("\n--- Time × Amb interaction ---\n")
if (length(ta_row) > 0) {
  e <- fixef_tab[ta_row, "Estimate"]
  s <- fixef_tab[ta_row, "Std. Error"]
  t <- fixef_tab[ta_row, "t value"]
  p <- fixef_tab[ta_row, "Pr(>|t|)"]
  cat(sprintf("  Est = %.4f, SE = %.4f, t = %.3f, p = %.4f\n", e, s, t, p))
} else {
  # fallback
  ta_row <- grep("time_c:amb_c", rownames(fixef_tab))
  ta_row <- ta_row[!grepl("cert", rownames(fixef_tab)[ta_row])]
  if (length(ta_row) > 0) {
    e <- fixef_tab[ta_row[1], "Estimate"]
    s <- fixef_tab[ta_row[1], "Std. Error"]
    t <- fixef_tab[ta_row[1], "t value"]
    p <- fixef_tab[ta_row[1], "Pr(>|t|)"]
    cat(sprintf("  Est = %.4f, SE = %.4f, t = %.3f, p = %.4f\n", e, s, t, p))
  }
}

# ——————————————————————————————————————————————————————————————
# D5. Variance components
# ——————————————————————————————————————————————————————————————

vc <- VarCorr(m_final)
cat(sprintf("\nICC = %.4f\n",
            as.numeric(vc$ID) / (as.numeric(vc$ID) + attr(vc, "sc")^2)))

# ——————————————————————————————————————————————————————————————
# D6. Conditional time slopes
# ——————————————————————————————————————————————————————————————

cat("\n--- Conditional time slopes at ±1 SD of cert & amb ---\n")

cert_sd <- sd(S2.full$cert, na.rm = TRUE)
amb_sd  <- sd(S2.full$amb,  na.rm = TRUE)

# 获取实际 time_c 值（中心化后的月份）
tc_vals <- sort(unique(long_data$time_c))

cond_grid <- expand.grid(
  time_c = tc_vals,  # T1, T2, T3 (centered months)
  cert_c = c(-cert_sd, 0, cert_sd),
  amb_c  = c(-amb_sd, 0, amb_sd),
  ex_c   = 0
)

cond_grid$pred <- predict(m_final, newdata = cond_grid, re.form = NA)

# Overall slope from T1 to T3 for each combo (per-month change × 12 months)
cond_slopes <- cond_grid %>%
  group_by(cert_c, amb_c) %>%
  summarise(
    change_T1_T3 = (pred[time_c == tc_vals[3]] - pred[time_c == tc_vals[1]]),
    .groups = "drop"
  ) %>%
  mutate(
    cert_lbl = case_when(
      cert_c == -cert_sd ~ "-1 SD", cert_c == 0 ~ "Mean", cert_c == cert_sd ~ "+1 SD"
    ),
    amb_lbl = case_when(
      amb_c == -amb_sd ~ "-1 SD", amb_c == 0 ~ "Mean", amb_c == amb_sd ~ "+1 SD"
    )
  )

print(as.data.frame(cond_slopes[, c("cert_lbl", "amb_lbl", "change_T1_T3")]))

cat("\n========== ALL ANALYSES COMPLETE ==========\n")
