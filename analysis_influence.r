##rm(list = ls(all = TRUE))
source('~/thib/projects/tools/R_lib.r')
setwd('/Users/thibault/thib/projects/reliable_info/reliable_info_influence/')

library(tidyverse)
library(posterior)
library(bayesplot)
library(loo)
library(knitr)

exp <- '6'
model_name <- 'linear_choice_influence'

# ---------------------------------------------------
## LOAD DATA AND FIT
# ---------------------------------------------------
load(paste0("./data/data_list_exp", exp, ".rdata"))
load(paste0("./data/data_reliability_exp", exp, ".rdata"))

data <- data %>%
  mutate(choice = case_when(
    (ResponseButtonOrder == 1 & Response == 0) ~ 2,
    (ResponseButtonOrder == 1 & Response == 1) ~ 1,
    (ResponseButtonOrder == 0 & Response == 0) ~ 1,
    (ResponseButtonOrder == 0 & Response == 1) ~ 2
  )) %>%
  mutate_at(vars(starts_with("color")), ~ ifelse(. == "blue", 1, 2)) %>%
  rowwise() %>%
  mutate(sample_number = sum(!is.na(c_across(starts_with("proba_"))))) %>%
  ungroup() %>%
  rename(influence_proba = SliderReliability, influence = SliderResponse) %>%
  mutate(influence = influence / 100 - 0.5)

data <- data %>%
  mutate(across(proba_1:proba_6, ~ as.integer(.x == influence_proba), .names = "{.col}_inf")) %>%
  rename_with(~ str_replace(.x, "proba_(\\d+)_inf", "influence_\\1"))

N <- data_list$N
T_max <- data_list$T_max
I_max <- data_list$I_max
subjs <- unique(data$ParticipantPrivateID)
t_subjs <- data_list$Tsubj

load(paste0("./results/fits/exp", exp, "/fit_", model_name, "_exp", exp, ".rdata"))

dir.create(paste0("./results/plots/exp", exp), recursive = TRUE, showWarnings = FALSE)
dir.create(paste0("./results/summary/exp", exp), recursive = TRUE, showWarnings = FALSE)

# ---------------------------------------------------
## EXTRACT ALL DRAWS ONCE
# ---------------------------------------------------
posterior_df <- as_draws_df(fit$draws())

# ---------------------------------------------------
## MODEL DIAGNOSIS
# ---------------------------------------------------

## Pairs plot
selected_params <- posterior_df[, grepl("^mu_", colnames(posterior_df)) & !grepl("^mu_pr", colnames(posterior_df))]
plot <- mcmc_pairs(selected_params)
ggsave(plot, file = paste0("./results/plots/exp", exp, "/pairs_plot_", model_name, "_exp", exp, ".pdf"))

## Trace plot
keep <- grepl('^mu_', colnames(posterior_df)) | grepl('^mu_pr\\[', colnames(posterior_df))
selected_params_trace <- posterior_df[, keep, drop = FALSE]
plot_trace <- mcmc_trace(selected_params_trace)
ggsave(filename = paste0("./results/plots/exp", exp, "/trace_plot_", model_name, "_exp", exp, ".pdf"),
       plot = plot_trace, width = 9, height = 6)

# ---------------------------------------------------
## MODEL SUMMARY
# ---------------------------------------------------
summary_stats <- posterior_summary(selected_params)
summary_df <- as.data.frame(summary_stats)

html_content <- kable(summary_df, format = "html", table.attr = "class='table table-bordered'")
latex_content <- kable(summary_df, format = "latex", booktabs = TRUE)

writeLines(html_content, paste0("./results/summary/exp", exp, "/summary_mu_", model_name, "_exp", exp, ".html"))
writeLines(latex_content, paste0("./results/summary/exp", exp, "/summary_mu_", model_name, "_exp", exp, ".tex"))

# ---------------------------------------------------
## EXTRACT PARAMETERS
# ---------------------------------------------------
invlogit <- function(x) 1 / (1 + exp(-x))
Phi_approx_R <- function(x) pnorm(x)

get_est <- function(sum_df, name) {
  if (!name %in% rownames(sum_df)) return(NA_real_)
  sum_df[name, "Estimate"]
}

mu_w1    <- get_est(summary_df, "mu_w1")
mu_w2    <- get_est(summary_df, "mu_w2")
mu_w3    <- get_est(summary_df, "mu_w3")
mu_w4    <- get_est(summary_df, "mu_w4")
mu_w5    <- get_est(summary_df, "mu_w5")
mu_alpha <- get_est(summary_df, "mu_alpha")
mu_beta  <- get_est(summary_df, "mu_beta")
mu_a     <- get_est(summary_df, "mu_a")
mu_b     <- get_est(summary_df, "mu_b")
mu_sigma_infl <- get_est(summary_df, "mu_sigma_infl")

cat("\n=== Group-level parameters ===\n")
cat(sprintf("mu_w1 = %.3f\nmu_w2 = %.3f\nmu_w3 = %.3f\nmu_w4 = %.3f\nmu_w5 = %.3f\n",
            mu_w1, mu_w2, mu_w3, mu_w4, mu_w5))
cat(sprintf("mu_alpha = %.3f\nmu_beta = %.3f\n", mu_alpha, mu_beta))
cat(sprintf("mu_a = %.3f\nmu_b = %.3f\nmu_sigma_infl = %.3f\n", mu_a, mu_b, mu_sigma_infl))

## Individual-level parameters
params_vars <- names(posterior_df)[grepl("^params\\[[0-9]+,[0-9]+\\]$", names(posterior_df))]
params_med <- vapply(params_vars, function(v) median(posterior_df[[v]], na.rm = TRUE), numeric(1))
idx <- stringr::str_match(params_vars, "^params\\[([0-9]+),([0-9]+)\\]$")

params_indiv <- tibble(
  subj = as.integer(idx[, 2]),
  k = as.integer(idx[, 3]),
  median = as.numeric(params_med)
) %>%
  mutate(k = paste0("k", k)) %>%
  pivot_wider(names_from = k, values_from = median) %>%
  arrange(subj) %>%
  rename(w1 = k1, w2 = k2, w3 = k3, w4 = k4, w5 = k5,
         alpha = k6, beta = k7, a_infl = k8, b_infl = k9, sigma_infl = k10)

# ---------------------------------------------------
## PROBABILITY TRANSFORMATION CURVES
# ---------------------------------------------------
p_grid <- (1:99) / 100
f_linear <- function(p, alpha, beta) invlogit(alpha * qlogis(p) + beta)

group_curve <- tibble(p = p_grid, fp = f_linear(p_grid, mu_alpha, mu_beta))

indiv_curves <- params_indiv %>%
  select(subj, alpha, beta) %>%
  crossing(p = p_grid) %>%
  mutate(fp = f_linear(p, alpha, beta))

plot_group <- ggplot() +
  geom_line(data = group_curve, aes(x = p, y = fp), linewidth = 1, color = "steelblue") +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
  coord_cartesian(xlim = c(0, 1), ylim = c(0, 1)) +
  theme_minimal() +
  labs(title = paste0("Probability transformation — group level (", model_name, ", exp", exp, ")"),
       x = "p", y = "f(p) = inv_logit(alpha * logit(p) + beta)")

ggsave(filename = paste0("./results/plots/exp", exp, "/proba_transform_group_", model_name, "_exp", exp, ".pdf"),
       plot = plot_group, width = 7, height = 5)

plot_indiv <- ggplot(indiv_curves, aes(x = p, y = fp, group = subj)) +
  geom_line(alpha = 0.15) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
  coord_cartesian(xlim = c(0, 1), ylim = c(0, 1)) +
  theme_minimal() +
  labs(title = paste0("Probability transformation — individuals (", model_name, ", exp", exp, ")"),
       x = "p", y = "f(p)")

ggsave(filename = paste0("./results/plots/exp", exp, "/proba_transform_indiv_spaghetti_", model_name, "_exp", exp, ".pdf"),
       plot = plot_indiv, width = 7, height = 5)

# ---------------------------------------------------
## SEQUENTIAL WEIGHTS PLOT
# ---------------------------------------------------
weights_group <- tibble(position = 1:6, weight = c(mu_w1, mu_w2, mu_w3, mu_w4, mu_w5, 1.0))

weights_indiv <- params_indiv %>%
  select(subj, w1, w2, w3, w4, w5) %>%
  mutate(w6 = 1.0) %>%
  pivot_longer(cols = w1:w6, names_to = "position", values_to = "weight") %>%
  mutate(position = as.integer(str_extract(position, "\\d+")))

plot_weights <- ggplot() +
  geom_line(data = weights_indiv, aes(x = position, y = weight, group = subj),
            alpha = 0.1, color = "grey50") +
  geom_line(data = weights_group, aes(x = position, y = weight),
            linewidth = 1.2, color = "steelblue") +
  geom_point(data = weights_group, aes(x = position, y = weight),
             size = 3, color = "steelblue") +
  geom_hline(yintercept = 1, linetype = "dashed", color = "red") +
  scale_x_continuous(breaks = 1:6, labels = paste0("w", 1:6)) +
  theme_minimal() +
  labs(title = paste0("Sequential weights (", model_name, ", exp", exp, ")"),
       x = "Sample position", y = "Weight")

ggsave(filename = paste0("./results/plots/exp", exp, "/sequential_weights_", model_name, "_exp", exp, ".pdf"),
       plot = plot_weights, width = 7, height = 5)

# ---------------------------------------------------
## POSTERIOR PREDICTIVE CHECK: CHOICE
# ---------------------------------------------------

pred_proba_vars <- names(posterior_df)[grepl("^pred_proba\\[", names(posterior_df))]
pred_proba_mean <- vapply(pred_proba_vars, function(v) mean(posterior_df[[v]], na.rm = TRUE), numeric(1))
idx_pp <- stringr::str_match(pred_proba_vars, "^pred_proba\\[([0-9]+),([0-9]+)\\]$")

df_pred_choice <- tibble(
  subj = as.integer(idx_pp[, 2]),
  trial = as.integer(idx_pp[, 3]),
  pred_p_blue = as.numeric(pred_proba_mean)
) %>% filter(pred_p_blue >= 0)

df_obs_choice <- data.frame()
for (n in 1:N) {
  for (t in 1:t_subjs[n]) {
    df_obs_choice <- rbind(df_obs_choice, data.frame(
      subj = n, trial = t, choice = data_list$choice[n, t]
    ))
  }
}

df_choice_check <- df_obs_choice %>%
  inner_join(df_pred_choice, by = c("subj", "trial")) %>%
  mutate(chose_blue = as.integer(choice == 1))

df_choice_check <- df_choice_check %>%
  mutate(pred_bin = cut(pred_p_blue, breaks = seq(0, 1, by = 0.1), include.lowest = TRUE))

df_calib <- df_choice_check %>%
  group_by(pred_bin) %>%
  summarise(
    obs_freq = mean(chose_blue),
    pred_mean = mean(pred_p_blue),
    n = n(),
    se = sqrt(obs_freq * (1 - obs_freq) / n),
    .groups = "drop"
  )

###############################################################################
## PLOT: Choice calibration
##
## DEFINITION: Predicted P(blue) is binned into deciles (0-10%, 10-20%, ...,
## 90-100%). Within each bin, x = mean predicted probability, y = observed
## frequency of blue choices, with 95% binomial CIs.
##
## INTERPRETATION: Points on the diagonal = well-calibrated model. Above the
## diagonal = underconfidence; below = overconfidence. Tests calibration only,
## not discrimination. A model predicting 50% everywhere would appear
## "calibrated" if the base rate is 50%.
###############################################################################
plot_calib <- ggplot(df_calib, aes(x = pred_mean, y = obs_freq)) +
  geom_pointrange(aes(ymin = obs_freq - 1.96 * se, ymax = obs_freq + 1.96 * se)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
  coord_cartesian(xlim = c(0, 1), ylim = c(0, 1)) +
  theme_minimal() +
  labs(title = paste0("Choice calibration (", model_name, ", exp", exp, ")"),
       x = "Predicted P(blue)", y = "Observed P(blue)")

ggsave(filename = paste0("./results/plots/exp", exp, "/choice_calibration_", model_name, "_exp", exp, ".pdf"),
       plot = plot_calib, width = 6, height = 5)

# ===================================================================
## POSTERIOR PREDICTIVE CHECK: INFLUENCE REPORT
# ===================================================================

## Extract predictions and delta_p
pred_infl_vars <- names(posterior_df)[grepl("^pred_influence\\[", names(posterior_df))]
pred_infl_mean <- vapply(pred_infl_vars, function(v) mean(posterior_df[[v]], na.rm = TRUE), numeric(1))
idx_pi <- stringr::str_match(pred_infl_vars, "^pred_influence\\[([0-9]+),([0-9]+)\\]$")

df_pred_infl <- tibble(
  subj = as.integer(idx_pi[, 2]),
  trial = as.integer(idx_pi[, 3]),
  pred_influence = as.numeric(pred_infl_mean)
) %>% filter(pred_influence > -90)

delta_p_vars <- names(posterior_df)[grepl("^delta_p_out\\[", names(posterior_df))]
delta_p_mean <- vapply(delta_p_vars, function(v) mean(posterior_df[[v]], na.rm = TRUE), numeric(1))
idx_dp <- stringr::str_match(delta_p_vars, "^delta_p_out\\[([0-9]+),([0-9]+)\\]$")

df_delta_p <- tibble(
  subj = as.integer(idx_dp[, 2]),
  trial = as.integer(idx_dp[, 3]),
  delta_p = as.numeric(delta_p_mean)
) %>% filter(delta_p > -90)

## Build observed influence data frame
df_obs_infl <- data.frame()
for (n in 1:N) {
  for (t in 1:t_subjs[n]) {
    infl_obs <- data_list$influence[n, t]
    if (infl_obs <= -90) next
    data_subj <- data %>% filter(ParticipantPrivateID == subjs[n])
    infl_proba <- data_subj$influence_proba[t]
    df_obs_infl <- rbind(df_obs_infl, data.frame(
      subj = n, trial = t, observed = infl_obs, influence_proba = infl_proba
    ))
  }
}

df_infl_check <- df_obs_infl %>%
  inner_join(df_pred_infl, by = c("subj", "trial")) %>%
  inner_join(df_delta_p, by = c("subj", "trial"))

## Add individual b_infl to the check dataframe for demeaning
df_infl_check <- df_infl_check %>%
  left_join(params_indiv %>% select(subj, a_infl, b_infl), by = "subj")

cat(sprintf("\n=== Influence check: %d rows after join ===\n", nrow(df_infl_check)))

# ---------------------------------------------------
## INFLUENCE PLOT 1: By reliability level
# ---------------------------------------------------
###############################################################################
## PLOT: Influence by reliability level
##
## DEFINITION: For each queried reliability level (e.g. 50%, 55%, 65%), the
## mean observed influence report (black) and mean predicted influence (red)
## are shown with 95% CIs, averaged across all subjects and trials.
##
## INTERPRETATION: Tests whether subjects report higher influence for more
## reliable samples and whether the model captures this pattern. This is
## the most direct test of the core model prediction. Note: predictions
## benefit from both a (sensitivity to delta_p) and b (intercept), so a
## good fit here does not by itself prove that delta_p is doing useful work.
###############################################################################
df_infl_by_proba <- df_infl_check %>%
  group_by(influence_proba) %>%
  summarise(
    obs_mean = mean(observed),
    obs_se = sd(observed) / sqrt(n()),
    pred_mean = mean(pred_influence),
    pred_se = sd(pred_influence) / sqrt(n()),
    n = n(),
    .groups = "drop"
  )

plot_by_proba <- ggplot(df_infl_by_proba, aes(x = factor(influence_proba))) +
  geom_pointrange(aes(y = obs_mean,
                      ymin = obs_mean - 1.96 * obs_se,
                      ymax = obs_mean + 1.96 * obs_se),
                  color = "black", position = position_nudge(x = -0.1)) +
  geom_pointrange(aes(y = pred_mean,
                      ymin = pred_mean - 1.96 * pred_se,
                      ymax = pred_mean + 1.96 * pred_se),
                  color = "red", position = position_nudge(x = 0.1)) +
  theme_minimal() +
  labs(title = paste0("Influence by reliability level: observed (black) vs predicted (red)\n(",
                      model_name, ", exp", exp, ")"),
       x = "Queried reliability level",
       y = "Influence report")

ggsave(filename = paste0("./results/plots/exp", exp, "/influence_by_proba_", model_name, "_exp", exp, ".pdf"),
       plot = plot_by_proba, width = 7, height = 5)

# ---------------------------------------------------
## INFLUENCE PLOT 2: Demeaned influence vs delta_p (pooled within-subject)
# ---------------------------------------------------
###############################################################################
## PLOT: Within-subject relationship between delta_p and influence report
##
## DEFINITION: For each trial, the observed influence report is demeaned
## within subject (subtract that subject's mean report across all trials).
## This removes all between-subject differences that the individual intercept
## b_infl can trivially absorb. The demeaned report is plotted against delta_p
## (also demeaned within subject for consistency).
##
## INTERPRETATION: This is the most diagnostic plot for assessing whether
## delta_p does real predictive work beyond what the individual intercept
## provides. A clear positive trend (or negative for below-chance
## reliabilities) means that within a given subject, trials where delta_p
## is larger are associated with higher influence reports. A flat cloud
## means delta_p adds no within-subject predictive power and all the
## apparent fit comes from the individual intercepts.
###############################################################################
df_demeaned <- df_infl_check %>%
  group_by(subj) %>%
  mutate(
    obs_demeaned = observed - mean(observed),
    delta_p_demeaned = delta_p - mean(delta_p)
  ) %>%
  ungroup()

plot_demeaned <- ggplot(df_demeaned, aes(x = delta_p_demeaned, y = obs_demeaned)) +
  geom_point(alpha = 0.1, size = 0.5) +
  geom_smooth(method = "lm", se = TRUE, color = "blue") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
  theme_minimal() +
  labs(title = paste0("Within-subject: demeaned influence vs demeaned delta_p\n(",
                      model_name, ", exp", exp, ")"),
       x = expression(Delta * p ~ "(demeaned within subject)"),
       y = "Influence report (demeaned within subject)")

ggsave(filename = paste0("./results/plots/exp", exp, "/influence_demeaned_vs_delta_p_", model_name, "_exp", exp, ".pdf"),
       plot = plot_demeaned, width = 7, height = 5)

# ---------------------------------------------------
## INFLUENCE PLOT 3: Within-subject correlations histogram
# ---------------------------------------------------
###############################################################################
## PLOT: Distribution of within-subject correlations between delta_p and
## observed influence report
##
## DEFINITION: For each subject, compute the Pearson correlation between
## delta_p and the observed influence report across their trials. Plot the
## distribution of these correlations as a histogram with the group mean
## marked.
##
## INTERPRETATION: If delta_p is a genuine within-subject predictor,
## most correlations should be positive (exp6) or negative (exp7).
## If correlations scatter around zero, delta_p is not tracking
## within-subject variability in reports. The spread indicates
## heterogeneity: some subjects may be more introspectively accurate
## than others.
###############################################################################
within_subj_cor <- df_infl_check %>%
  group_by(subj) %>%
  summarise(
    r = cor(delta_p, observed, use = "complete.obs"),
    n_trials = n(),
    .groups = "drop"
  ) %>%
  filter(!is.na(r))

mean_r <- mean(within_subj_cor$r)

plot_cor_hist <- ggplot(within_subj_cor, aes(x = r)) +
  geom_histogram(bins = 20, fill = "steelblue", alpha = 0.7, color = "white") +
  geom_vline(xintercept = mean_r, color = "red", linewidth = 1, linetype = "dashed") +
  geom_vline(xintercept = 0, color = "grey50", linetype = "dashed") +
  annotate("text", x = mean_r, y = Inf, label = sprintf("mean r = %.3f", mean_r),
           vjust = 2, hjust = -0.1, color = "red", size = 4) +
  theme_minimal() +
  labs(title = paste0("Within-subject correlations: delta_p vs influence report\n(",
                      model_name, ", exp", exp, ")"),
       x = "Pearson r (delta_p, influence report)",
       y = "Number of subjects")

ggsave(filename = paste0("./results/plots/exp", exp, "/within_subject_cor_hist_", model_name, "_exp", exp, ".pdf"),
       plot = plot_cor_hist, width = 7, height = 5)

# ---------------------------------------------------
## INFLUENCE PLOT 4: Estimated a_infl vs empirical within-subject slope
# ---------------------------------------------------
###############################################################################
## PLOT: Model-estimated a_infl vs empirical within-subject regression slope
##
## DEFINITION: For each subject, fit a simple linear regression of observed
## influence on delta_p (empirical slope). Compare this to the model's
## estimated a_infl parameter. If the model correctly captures individual
## differences in sensitivity to delta_p, these should correlate.
##
## INTERPRETATION: Points along the diagonal = the hierarchical model
## recovers individual sensitivities well. A positive correlation with
## slope < 1 indicates hierarchical shrinkage (expected). A flat cloud
## means the model's a_infl does not track the empirical relationship,
## suggesting misspecification or that delta_p is too noisy to estimate
## individual slopes reliably.
###############################################################################
empirical_slopes <- df_infl_check %>%
  group_by(subj) %>%
  summarise(
    emp_slope = coef(lm(observed ~ delta_p))[2],
    .groups = "drop"
  )

slopes_vs_a <- empirical_slopes %>%
  inner_join(params_indiv %>% select(subj, a_infl), by = "subj")

plot_slopes <- ggplot(slopes_vs_a, aes(x = emp_slope, y = a_infl)) +
  geom_point(size = 2) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
  geom_smooth(method = "lm", se = TRUE, color = "blue") +
  theme_minimal() +
  labs(title = paste0("Empirical slope vs estimated a_infl (", model_name, ", exp", exp, ")"),
       x = "Empirical within-subject slope (influence ~ delta_p)",
       y = "Model-estimated a_infl")

ggsave(filename = paste0("./results/plots/exp", exp, "/empirical_slope_vs_a_infl_", model_name, "_exp", exp, ".pdf"),
       plot = plot_slopes, width = 6, height = 5)

# ---------------------------------------------------
## INFLUENCE PLOT 5: Demeaned density (within-subject variability)
# ---------------------------------------------------
###############################################################################
## PLOT: Within-subject variability: observed vs predicted (demeaned)
##
## DEFINITION: For each subject, compute their mean observed influence and
## their mean predicted influence (from one posterior draw) across all trials.
## Subtract these means from trial-level values, yielding "deviations from
## subject mean" for both observed and predicted. Overlay the two densities.
##
## This is model-free demeaning: it removes all between-subject differences
## regardless of whether they are captured by b_infl, by a_infl times the
## subject's average delta_p, or by anything else. What remains is purely
## the within-subject trial-to-trial variability.
##
## The predicted distribution uses a single posterior draw (not the posterior
## mean) to preserve the full noise structure from sigma_infl. This makes
## the comparison fair: both distributions include trial-level noise.
##
## INTERPRETATION: If the two distributions overlap, the model captures
## within-subject variability well — meaning a_infl * delta_p plus noise
## reproduces the trial-to-trial fluctuations. If the predicted distribution
## is too narrow, the model underpredicts within-subject variability (sigma
## too small or delta_p has insufficient range). If too wide, the model
## overpredicts variability.
###############################################################################

## Get one full posterior draw of pred_influence
pred_infl_mat <- as.matrix(posterior_df[, pred_infl_vars])
one_draw_idx <- sample(1:nrow(pred_infl_mat), 1)
one_draw <- pred_infl_mat[one_draw_idx, ]

df_one_draw <- tibble(
  subj = as.integer(idx_pi[, 2]),
  trial = as.integer(idx_pi[, 3]),
  pred_draw = as.numeric(one_draw)
) %>% filter(pred_draw > -90)

## Join with observed and demean within subject
df_demeaned_dens <- df_infl_check %>%
  select(subj, trial, observed) %>%
  inner_join(df_one_draw, by = c("subj", "trial")) %>%
  group_by(subj) %>%
  mutate(
    obs_demeaned = observed - mean(observed),
    pred_demeaned = pred_draw - mean(pred_draw)
  ) %>%
  ungroup()

df_dens_both <- bind_rows(
  df_demeaned_dens %>% select(value = obs_demeaned) %>% mutate(type = "Observed (demeaned)"),
  df_demeaned_dens %>% select(value = pred_demeaned) %>% mutate(type = "Predicted (demeaned, 1 draw)")
)

plot_demeaned_dens <- ggplot(df_dens_both, aes(x = value, fill = type)) +
  geom_density(alpha = 0.4) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
  theme_minimal() +
  labs(title = paste0("Within-subject variability: observed vs predicted (demeaned)\n(",
                      model_name, ", exp", exp, ")"),
       x = "Deviation from subject mean", y = "Density",
       fill = "")

ggsave(filename = paste0("./results/plots/exp", exp, "/demeaned_density_", model_name, "_exp", exp, ".pdf"),
       plot = plot_demeaned_dens, width = 7, height = 5)

# ---------------------------------------------------
## INFLUENCE PLOT 6: delta_p vs observed, faceted by reliability
# ---------------------------------------------------
###############################################################################
## PLOT: delta_p vs observed influence, faceted by queried reliability level
##
## DEFINITION: Same as the demeaned scatter but NOT demeaned, split by
## reliability level. Shows the raw relationship between the model's
## delta_p predictor and observed reports for each reliability level.
##
## INTERPRETATION: Checks whether the delta_p -> influence relationship
## holds consistently across reliability levels. Different slopes across
## facets suggest subjects use different introspective strategies for
## different reliability levels. Also reveals whether delta_p has
## sufficient range at each level (narrow range = little signal).
###############################################################################
plot_delta_p_facet <- ggplot(df_infl_check, aes(x = delta_p, y = observed)) +
  geom_point(alpha = 0.15, size = 0.5) +
  geom_smooth(method = "lm", se = TRUE, color = "blue") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
  facet_wrap(~ influence_proba) +
  theme_minimal() +
  labs(title = paste0("Delta P vs observed report by reliability level\n(",
                      model_name, ", exp", exp, ")"),
       x = expression(Delta * p),
       y = "Observed influence report")

ggsave(filename = paste0("./results/plots/exp", exp, "/delta_p_vs_observed_by_proba_", model_name, "_exp", exp, ".pdf"),
       plot = plot_delta_p_facet, width = 9, height = 6)

# ---------------------------------------------------
## INFLUENCE PLOT 7: Influence parameters (a, b) across subjects
# ---------------------------------------------------
###############################################################################
## PLOT: Individual influence parameters a_infl vs b_infl
##
## DEFINITION: One point per subject in the (a_infl, b_infl) space. Red
## dashed lines mark group-level means.
##
## INTERPRETATION: a_infl = sensitivity of report to delta_p. Positive a
## means reports track the counterfactual shift. b_infl = baseline reporting
## tendency. The spread reveals heterogeneity in introspective ability and
## reporting bias. Correlations between a and b may suggest trade-offs or
## reparameterisation opportunities.
###############################################################################
plot_ab <- ggplot(params_indiv, aes(x = a_infl, y = b_infl)) +
  geom_point(alpha = 0.5) +
  geom_vline(xintercept = mu_a, color = "red", linetype = "dashed") +
  geom_hline(yintercept = mu_b, color = "red", linetype = "dashed") +
  theme_minimal() +
  labs(title = paste0("Influence parameters across subjects (", model_name, ", exp", exp, ")"),
       x = "a (scaling)", y = "b (shift)")

ggsave(filename = paste0("./results/plots/exp", exp, "/influence_params_ab_", model_name, "_exp", exp, ".pdf"),
       plot = plot_ab, width = 6, height = 5)

# ---------------------------------------------------
## INFLUENCE PLOT 8: Delta_p density
# ---------------------------------------------------
###############################################################################
## PLOT: Distribution of delta_p values
##
## DEFINITION: Density of the model-computed delta_p across all trials and
## subjects. delta_p = P(choice|all) - P(choice|without queried samples).
##
## INTERPRETATION: Shows the range and shape of the key influence predictor.
## Concentrated near zero = queried samples rarely mattered (remaining
## evidence was decisive). Wide spread = good variability for the influence
## model to exploit. If concentrated, the influence model has little signal
## regardless of parameter quality.
###############################################################################
plot_delta_p_dens <- ggplot(df_infl_check, aes(x = delta_p)) +
  geom_density(fill = "steelblue", alpha = 0.4) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
  theme_minimal() +
  labs(title = paste0("Distribution of delta P (", model_name, ", exp", exp, ")"),
       x = expression(Delta * p),
       y = "Density")

ggsave(filename = paste0("./results/plots/exp", exp, "/delta_p_density_", model_name, "_exp", exp, ".pdf"),
       plot = plot_delta_p_dens, width = 7, height = 5)

cat("\n=== Done. All plots saved. ===\n")
