###############################################################################
##  analysis_influence_tools.r
##
##  Helper functions for analysis_influence.r.
##  Source this file at the top of the main script.
##
##  Contents:
##    1. Utilities
##    2. Posterior loading (handles both original and renamed fits)
##    3. Parameter extraction
##    4. Data-frame builders
##    5. Plot functions (diagnostics, choice, influence, LOO)
###############################################################################

library(tidyverse)
library(posterior)
library(bayesplot)
library(loo)
library(knitr)

# ===================================================================
# 1. UTILITIES
# ===================================================================

invlogit <- function(x) 1 / (1 + exp(-x))

## Safely extract a value from a summary data.frame by row name
get_est <- function(sum_df, name) {
  if (!name %in% rownames(sum_df)) return(NA_real_)
  sum_df[name, "Estimate"]
}

# ===================================================================
# 2. POSTERIOR LOADING
# ===================================================================

#' Load posterior draws for a model.
#'
#' If rename_lo_fit.r has been run, the lo .rdata file contains
#' 'lo_draws' (a draws_df with renamed variables). We use that
#' when available; otherwise we fall back to fit$draws().
load_posterior <- function(model_name, exp) {
  fit_path <- paste0("./results/fits/exp", exp, "/fit_", model_name, "_exp", exp, ".rdata")
  env <- new.env()
  load(fit_path, envir = env)

  ## Otherwise extract from the fit object
  cat(sprintf("  Extracting draws from fit for %s\n", model_name))
  as_draws_df(env$fit$draws())
}

# ===================================================================
# 3. PARAMETER EXTRACTION
# ===================================================================

#' Extract group-level means from posterior draws.
#'
#' After renaming, both models have mu_w1..mu_w5, mu_alpha, mu_beta,
#' mu_a_infl, mu_b_infl, mu_sigma_infl as named variables.
extract_group_means <- function(posterior_df) {
  get_mean <- function(name) {
    if (name %in% names(posterior_df)) mean(posterior_df[[name]], na.rm = TRUE)
    else NA_real_
  }
  list(
    mu_w1    = get_mean("mu_w1"),
    mu_w2    = get_mean("mu_w2"),
    mu_w3    = get_mean("mu_w3"),
    mu_w4    = get_mean("mu_w4"),
    mu_w5    = get_mean("mu_w5"),
    mu_alpha = get_mean("mu_alpha"),
    mu_beta  = get_mean("mu_beta"),
    mu_a     = get_mean("mu_a_infl"),
    mu_b     = get_mean("mu_b_infl"),
    mu_sigma = get_mean("mu_sigma_infl")
  )
}

#' Extract individual-level parameters (median across draws).
#'
#' Both models store subject params in params[n, k] with k = 1..10.
extract_indiv_params <- function(posterior_df) {
  params_vars <- names(posterior_df)[grepl("^params\\[[0-9]+,[0-9]+\\]$", names(posterior_df))]
  params_med  <- vapply(params_vars, function(v) median(posterior_df[[v]], na.rm = TRUE), numeric(1))
  idx <- stringr::str_match(params_vars, "^params\\[([0-9]+),([0-9]+)\\]$")

  tibble(
    subj   = as.integer(idx[, 2]),
    k      = as.integer(idx[, 3]),
    median = as.numeric(params_med)
  ) %>%
    mutate(k = paste0("k", k)) %>%
    pivot_wider(names_from = k, values_from = median) %>%
    arrange(subj) %>%
    rename(w1 = k1, w2 = k2, w3 = k3, w4 = k4, w5 = k5,
           alpha = k6, beta = k7, a_infl = k8, b_infl = k9, sigma_infl = k10)
}

# ===================================================================
# 4. DATA-FRAME BUILDERS
# ===================================================================

#' Build a data.frame of observed choices.
build_obs_choice_df <- function(data_list, t_subjs, N) {
  rows <- vector("list", sum(t_subjs))
  i <- 1
  for (n in 1:N) {
    for (t in 1:t_subjs[n]) {
      rows[[i]] <- data.frame(subj = n, trial = t, choice = data_list$choice[n, t])
      i <- i + 1
    }
  }
  bind_rows(rows)
}

#' Build a data.frame of observed influence reports.
#' Skips trials with missing influence (< -90).
build_obs_influence_df <- function(data_list, data, subjs, t_subjs, N) {
  rows <- list()
  for (n in 1:N) {
    data_subj <- data %>% filter(ParticipantPrivateID == subjs[n])
    for (t in 1:t_subjs[n]) {
      infl_obs <- data_list$influence[n, t]
      if (infl_obs <= -90) next
      rows[[length(rows) + 1]] <- data.frame(
        subj = n, trial = t,
        observed = infl_obs,
        influence_proba = data_subj$influence_proba[t]
      )
    }
  }
  bind_rows(rows)
}

#' Extract posterior-mean predicted P(blue) per trial.
extract_pred_choice <- function(posterior_df) {
  vars <- names(posterior_df)[grepl("^pred_proba\\[", names(posterior_df))]
  vals <- vapply(vars, function(v) mean(posterior_df[[v]], na.rm = TRUE), numeric(1))
  idx  <- stringr::str_match(vars, "^pred_proba\\[([0-9]+),([0-9]+)\\]$")
  tibble(
    subj        = as.integer(idx[, 2]),
    trial       = as.integer(idx[, 3]),
    pred_p_blue = as.numeric(vals)
  ) %>% filter(pred_p_blue >= 0)
}

#' Extract posterior-mean predicted influence per trial.
extract_pred_influence <- function(posterior_df) {
  vars <- names(posterior_df)[grepl("^pred_influence\\[", names(posterior_df))]
  vals <- vapply(vars, function(v) mean(posterior_df[[v]], na.rm = TRUE), numeric(1))
  idx  <- stringr::str_match(vars, "^pred_influence\\[([0-9]+),([0-9]+)\\]$")
  tibble(
    subj           = as.integer(idx[, 2]),
    trial          = as.integer(idx[, 3]),
    pred_influence = as.numeric(vals)
  ) %>% filter(pred_influence > -90)
}

#' Extract the influence predictor (sub_ev_out or delta_p_out).
extract_influence_predictor <- function(posterior_df, var_prefix) {
  pattern <- paste0("^", var_prefix, "\\[")
  vars <- names(posterior_df)[grepl(pattern, names(posterior_df))]
  if (length(vars) == 0) {
    warning(sprintf("No variables matching '%s' found in posterior draws", var_prefix))
    return(tibble(subj = integer(), trial = integer(), predictor = numeric()))
  }
  vals <- vapply(vars, function(v) mean(posterior_df[[v]], na.rm = TRUE), numeric(1))
  idx  <- stringr::str_match(vars, paste0("^", var_prefix, "\\[([0-9]+),([0-9]+)\\]$"))
  tibble(
    subj      = as.integer(idx[, 2]),
    trial     = as.integer(idx[, 3]),
    predictor = as.numeric(vals)
  ) %>% filter(predictor > -90)
}

#' Build the full influence-check data.frame.
#' Joins observed reports, predicted influence, predictor values,
#' and individual parameters.
build_infl_check_df <- function(posterior_df, var_prefix,
                                data_list, data, subjs, t_subjs, N,
                                params_indiv) {
  df_obs  <- build_obs_influence_df(data_list, data, subjs, t_subjs, N)
  df_pred <- extract_pred_influence(posterior_df)
  df_iv   <- extract_influence_predictor(posterior_df, var_prefix)

  df_obs %>%
    inner_join(df_pred, by = c("subj", "trial")) %>%
    inner_join(df_iv,   by = c("subj", "trial")) %>%
    left_join(params_indiv %>% select(subj, a_infl, b_infl), by = "subj")
}

# ===================================================================
# 5. PLOT FUNCTIONS
# ===================================================================

## ---------- Diagnostics ----------

plot_pairs <- function(posterior_df, model_name, exp, out_dir) {
  mu_vars <- c("mu_w1","mu_w2","mu_w3","mu_w4","mu_w5",
               "mu_alpha","mu_beta","mu_a_infl","mu_b_infl","mu_sigma_infl")
  sel <- posterior_df[, intersect(mu_vars, colnames(posterior_df)), drop = FALSE]
  p <- mcmc_pairs(sel)
  ggsave(p, file = file.path(out_dir, paste0("pairs_plot_", model_name, "_exp", exp, ".pdf")))
}

plot_traces <- function(posterior_df, model_name, exp, out_dir) {
  keep <- grepl('^mu_pr\\[', colnames(posterior_df))
  sel  <- posterior_df[, keep, drop = FALSE]
  p    <- mcmc_trace(sel)
  ggsave(filename = file.path(out_dir, paste0("trace_plot_", model_name, "_exp", exp, ".pdf")),
         plot = p, width = 9, height = 6)
}

## ---------- Model summary ----------

save_summary <- function(posterior_df, model_name, exp, sum_dir) {
  mu_vars <- c("mu_w1","mu_w2","mu_w3","mu_w4","mu_w5",
               "mu_alpha","mu_beta","mu_a_infl","mu_b_infl","mu_sigma_infl")
  sel <- posterior_df[, intersect(mu_vars, colnames(posterior_df)), drop = FALSE]
  summary_tbl <- as.data.frame(posterior_summary(sel))

  writeLines(kable(summary_tbl, format = "html", table.attr = "class='table table-bordered'"),
             file.path(sum_dir, paste0("summary_", model_name, "_exp", exp, ".html")))
  writeLines(kable(summary_tbl, format = "latex", booktabs = TRUE),
             file.path(sum_dir, paste0("summary_", model_name, "_exp", exp, ".tex")))
  summary_tbl
}

## ---------- Choice ----------

plot_proba_transform <- function(mu_alpha, mu_beta, params_indiv,
                                  model_name, exp, out_dir) {
  p_grid <- (1:99) / 100
  f_lin  <- function(p, a, b) invlogit(a * qlogis(p) + b)

  group_curve  <- tibble(p = p_grid, fp = f_lin(p_grid, mu_alpha, mu_beta))
  indiv_curves <- params_indiv %>%
    select(subj, alpha, beta) %>%
    crossing(p = p_grid) %>%
    mutate(fp = f_lin(p, alpha, beta))

  p1 <- ggplot() +
    geom_line(data = group_curve, aes(x = p, y = fp), linewidth = 1, color = "steelblue") +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
    coord_cartesian(xlim = c(0, 1), ylim = c(0, 1)) +
    theme_minimal() +
    labs(title = paste0("Probability transformation — group (", model_name, ", exp", exp, ")"),
         x = "p", y = "f(p)")
  ggsave(file.path(out_dir, paste0("proba_transform_group_", model_name, "_exp", exp, ".pdf")),
         plot = p1, width = 7, height = 5)

  p2 <- ggplot(indiv_curves, aes(x = p, y = fp, group = subj)) +
    geom_line(alpha = 0.15) +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
    coord_cartesian(xlim = c(0, 1), ylim = c(0, 1)) +
    theme_minimal() +
    labs(title = paste0("Probability transformation — individuals (", model_name, ", exp", exp, ")"),
         x = "p", y = "f(p)")
  ggsave(file.path(out_dir, paste0("proba_transform_indiv_", model_name, "_exp", exp, ".pdf")),
         plot = p2, width = 7, height = 5)
}

plot_seq_weights <- function(gm, params_indiv, model_name, exp, out_dir) {
  weights_group <- tibble(position = 1:6,
                          weight = c(gm$mu_w1, gm$mu_w2, gm$mu_w3,
                                     gm$mu_w4, gm$mu_w5, 1.0))
  weights_indiv <- params_indiv %>%
    select(subj, w1, w2, w3, w4, w5) %>%
    mutate(w6 = 1.0) %>%
    pivot_longer(cols = w1:w6, names_to = "position", values_to = "weight") %>%
    mutate(position = as.integer(str_extract(position, "\\d+")))

  p <- ggplot() +
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
  ggsave(file.path(out_dir, paste0("sequential_weights_", model_name, "_exp", exp, ".pdf")),
         plot = p, width = 7, height = 5)
}

plot_choice_calibration <- function(posterior_df, data_list, t_subjs, N,
                                    model_name, exp, out_dir) {
  df_pred <- extract_pred_choice(posterior_df)
  df_obs  <- build_obs_choice_df(data_list, t_subjs, N)

  df <- df_obs %>%
    inner_join(df_pred, by = c("subj", "trial")) %>%
    mutate(chose_blue = as.integer(choice == 1),
           pred_bin   = cut(pred_p_blue, breaks = seq(0, 1, by = 0.1), include.lowest = TRUE))

  df_calib <- df %>%
    group_by(pred_bin) %>%
    summarise(obs_freq  = mean(chose_blue),
              pred_mean = mean(pred_p_blue),
              n = n(),
              se = sqrt(obs_freq * (1 - obs_freq) / n),
              .groups = "drop")

  p <- ggplot(df_calib, aes(x = pred_mean, y = obs_freq)) +
    geom_pointrange(aes(ymin = obs_freq - 1.96 * se, ymax = obs_freq + 1.96 * se)) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
    coord_cartesian(xlim = c(0, 1), ylim = c(0, 1)) +
    theme_minimal() +
    labs(title = paste0("Choice calibration (", model_name, ", exp", exp, ")"),
         x = "Predicted P(blue)", y = "Observed P(blue)")
  ggsave(file.path(out_dir, paste0("choice_calibration_", model_name, "_exp", exp, ".pdf")),
         plot = p, width = 6, height = 5)
}

#' Per-subject choice calibration (one panel per subject).
#'
#' Each facet shows one subject's calibration: predicted P(blue) binned
#' into quintiles vs observed frequency of choosing blue.
plot_choice_calibration_per_subj <- function(posterior_df, data_list, t_subjs, N,
                                              model_name, exp, out_dir) {
  df_pred <- extract_pred_choice(posterior_df)
  df_obs  <- build_obs_choice_df(data_list, t_subjs, N)

  df <- df_obs %>%
    inner_join(df_pred, by = c("subj", "trial")) %>%
    mutate(chose_blue = as.integer(choice == 1))

  ## Bin within each subject (quintiles)
  df <- df %>%
    group_by(subj) %>%
    mutate(pred_bin = cut(pred_p_blue,
                          breaks = quantile(pred_p_blue, probs = seq(0, 1, by = 0.2)),
                          include.lowest = TRUE)) %>%
    ungroup() %>%
    filter(!is.na(pred_bin))

  df_calib <- df %>%
    group_by(subj, pred_bin) %>%
    summarise(obs_freq  = mean(chose_blue),
              pred_mean = mean(pred_p_blue),
              n = n(),
              se = sqrt(obs_freq * (1 - obs_freq) / n),
              .groups = "drop")

  p <- ggplot(df_calib, aes(x = pred_mean, y = obs_freq)) +
    geom_pointrange(aes(ymin = pmax(obs_freq - 1.96 * se, 0),
                        ymax = pmin(obs_freq + 1.96 * se, 1)),
                    size = 0.3) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
    coord_cartesian(xlim = c(0, 1), ylim = c(0, 1)) +
    facet_wrap(~ subj) +
    theme_minimal(base_size = 8) +
    labs(title = paste0("Choice calibration per subject (", model_name, ", exp", exp, ")"),
         x = "Predicted P(blue)", y = "Observed P(blue)")
  ggsave(file.path(out_dir, paste0("choice_calibration_per_subj_", model_name, "_exp", exp, ".pdf")),
         plot = p, width = 14, height = 10)
}

## ---------- Influence plots ----------

plot_infl_by_reliability <- function(df, model_name, exp, out_dir) {
  df_agg <- df %>%
    group_by(influence_proba) %>%
    summarise(obs_mean = mean(observed), obs_se = sd(observed) / sqrt(n()),
              pred_mean = mean(pred_influence), pred_se = sd(pred_influence) / sqrt(n()),
              .groups = "drop") %>%
    mutate(x_num = as.numeric(factor(influence_proba)))

  p <- ggplot(df_agg) +
    geom_line(aes(x = x_num - 0.1, y = obs_mean), color = "black") +
    geom_pointrange(aes(x = x_num - 0.1, y = obs_mean,
                        ymin = obs_mean - 1.96*obs_se, ymax = obs_mean + 1.96*obs_se),
                    color = "black") +
    geom_line(aes(x = x_num + 0.1, y = pred_mean), color = "red") +
    geom_pointrange(aes(x = x_num + 0.1, y = pred_mean,
                        ymin = pred_mean - 1.96*pred_se, ymax = pred_mean + 1.96*pred_se),
                    color = "red") +
    scale_x_continuous(breaks = df_agg$x_num, labels = df_agg$influence_proba) +
    theme_minimal() +
    labs(title = paste0("Influence by reliability (", model_name, ", exp", exp, ")"),
         x = "Reliability level", y = "Influence report")
  ggsave(file.path(out_dir, paste0("influence_by_proba_", model_name, "_exp", exp, ".pdf")),
         plot = p, width = 7, height = 5)
}

#' Per-subject influence by reliability (one panel per subject).
#'
#' Each facet shows one subject's observed (black) vs predicted (red)
#' mean influence at each reliability level.
plot_infl_by_reliability_per_subj <- function(df, model_name, exp, out_dir) {

  df_subj <- df %>%
    group_by(subj, influence_proba) %>%
    summarise(obs_mean  = mean(observed),
              obs_se    = sd(observed) / sqrt(n()),
              pred_mean = mean(pred_influence),
              pred_se   = sd(pred_influence) / sqrt(n()),
              .groups = "drop") %>%
    mutate(x_num = as.numeric(factor(influence_proba)))

  x_labels <- df_subj %>% distinct(influence_proba, x_num) %>% arrange(x_num)

  p <- ggplot(df_subj) +
    geom_line(aes(x = x_num - 0.15, y = obs_mean), color = "black", linewidth = 0.3) +
    geom_pointrange(aes(x = x_num - 0.15, y = obs_mean,
                        ymin = obs_mean - 1.96 * obs_se,
                        ymax = obs_mean + 1.96 * obs_se),
                    color = "black", size = 0.2) +
    geom_line(aes(x = x_num + 0.15, y = pred_mean), color = "red", linewidth = 0.3) +
    geom_pointrange(aes(x = x_num + 0.15, y = pred_mean,
                        ymin = pred_mean - 1.96 * pred_se,
                        ymax = pred_mean + 1.96 * pred_se),
                    color = "red", size = 0.2) +
    scale_x_continuous(breaks = x_labels$x_num, labels = x_labels$influence_proba) +
    facet_wrap(~ subj) +
    theme_minimal(base_size = 8) +
    labs(title = paste0("Influence by reliability per subject: obs (black) vs pred (red)\n(",
                        model_name, ", exp", exp, ")"),
         x = "Reliability level", y = "Influence report")
  ggsave(file.path(out_dir, paste0("influence_by_proba_per_subj_", model_name, "_exp", exp, ".pdf")),
         plot = p, width = 14, height = 10)
}

plot_infl_demeaned_scatter <- function(df, pred_label, model_name, exp, out_dir) {
  df_dm <- df %>%
    group_by(subj) %>%
    mutate(obs_dm = observed - mean(observed), pred_dm = predictor - mean(predictor)) %>%
    ungroup()

  p <- ggplot(df_dm, aes(x = pred_dm, y = obs_dm)) +
    geom_point(alpha = 0.1, size = 0.5) +
    geom_smooth(method = "lm", se = TRUE, color = "blue") +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
    theme_minimal() +
    labs(title = paste0("Demeaned influence vs ", pred_label, " (", model_name, ", exp", exp, ")"),
         x = paste0(pred_label, " (demeaned)"), y = "Influence (demeaned)")
  ggsave(file.path(out_dir, paste0("influence_demeaned_scatter_", model_name, "_exp", exp, ".pdf")),
         plot = p, width = 7, height = 5)
}

plot_infl_within_cor <- function(df, pred_label, model_name, exp, out_dir) {
  cors <- df %>%
    group_by(subj) %>%
    summarise(r = cor(predictor, observed, use = "complete.obs"), .groups = "drop") %>%
    filter(!is.na(r))
  mean_r <- mean(cors$r)

  p <- ggplot(cors, aes(x = r)) +
    geom_histogram(bins = 20, fill = "steelblue", alpha = 0.7, color = "white") +
    geom_vline(xintercept = mean_r, color = "red", linewidth = 1, linetype = "dashed") +
    geom_vline(xintercept = 0, color = "grey50", linetype = "dashed") +
    annotate("text", x = mean_r, y = Inf, label = sprintf("mean r = %.3f", mean_r),
             vjust = 2, hjust = -0.1, color = "red", size = 4) +
    theme_minimal() +
    labs(title = paste0("Within-subj cor: ", pred_label, " vs influence (", model_name, ", exp", exp, ")"),
         x = paste0("Pearson r"), y = "N subjects")
  ggsave(file.path(out_dir, paste0("within_subject_cor_", model_name, "_exp", exp, ".pdf")),
         plot = p, width = 7, height = 5)
}

plot_infl_slope_recovery <- function(df, pred_label, model_name, exp, out_dir) {
  emp <- df %>%
    group_by(subj) %>%
    summarise(emp_slope = coef(lm(observed ~ predictor))[2], .groups = "drop")
  merged <- emp %>% inner_join(df %>% distinct(subj, a_infl), by = "subj")

  p <- ggplot(merged, aes(x = emp_slope, y = a_infl)) +
    geom_point(size = 2) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
    geom_smooth(method = "lm", se = TRUE, color = "blue") +
    theme_minimal() +
    labs(title = paste0("Empirical slope vs a_infl (", model_name, ", exp", exp, ")"),
         x = paste0("Empirical slope (influence ~ ", pred_label, ")"), y = "Model a_infl")
  ggsave(file.path(out_dir, paste0("slope_recovery_", model_name, "_exp", exp, ".pdf")),
         plot = p, width = 6, height = 5)
}

plot_infl_demeaned_density <- function(df, posterior_df, model_name, exp, out_dir) {
  pred_infl_vars <- names(posterior_df)[grepl("^pred_influence\\[", names(posterior_df))]
  pred_infl_mat  <- as.matrix(posterior_df[, pred_infl_vars])
  one_draw       <- pred_infl_mat[sample(1:nrow(pred_infl_mat), 1), ]
  idx_pi         <- stringr::str_match(pred_infl_vars, "^pred_influence\\[([0-9]+),([0-9]+)\\]$")

  df_draw <- tibble(subj = as.integer(idx_pi[, 2]), trial = as.integer(idx_pi[, 3]),
                    pred_draw = as.numeric(one_draw)) %>%
    filter(pred_draw > -90)

  df_dm <- df %>% select(subj, trial, observed) %>%
    inner_join(df_draw, by = c("subj", "trial")) %>%
    group_by(subj) %>%
    mutate(obs_dm = observed - mean(observed), pred_dm = pred_draw - mean(pred_draw)) %>%
    ungroup()

  df_both <- bind_rows(
    df_dm %>% transmute(value = obs_dm,  type = "Observed (demeaned)"),
    df_dm %>% transmute(value = pred_dm, type = "Predicted (demeaned, 1 draw)")
  )

  p <- ggplot(df_both, aes(x = value, fill = type)) +
    geom_density(alpha = 0.4) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
    theme_minimal() +
    labs(title = paste0("Within-subj variability (", model_name, ", exp", exp, ")"),
         x = "Deviation from subject mean", y = "Density", fill = "")
  ggsave(file.path(out_dir, paste0("demeaned_density_", model_name, "_exp", exp, ".pdf")),
         plot = p, width = 7, height = 5)
}

plot_infl_by_proba_facet <- function(df, pred_label, model_name, exp, out_dir) {
  p <- ggplot(df, aes(x = predictor, y = observed)) +
    geom_point(alpha = 0.15, size = 0.5) +
    geom_smooth(method = "lm", se = TRUE, color = "blue") +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
    facet_wrap(~ influence_proba) +
    theme_minimal() +
    labs(title = paste0(pred_label, " vs observed by reliability (", model_name, ", exp", exp, ")"),
         x = pred_label, y = "Observed influence")
  ggsave(file.path(out_dir, paste0("predictor_vs_obs_facet_", model_name, "_exp", exp, ".pdf")),
         plot = p, width = 9, height = 6)
}

plot_infl_ab_scatter <- function(gm, params_indiv, model_name, exp, out_dir) {
  p <- ggplot(params_indiv, aes(x = a_infl, y = b_infl)) +
    geom_point(alpha = 0.5) +
    geom_vline(xintercept = gm$mu_a, color = "red", linetype = "dashed") +
    geom_hline(yintercept = gm$mu_b, color = "red", linetype = "dashed") +
    theme_minimal() +
    labs(title = paste0("Influence params (", model_name, ", exp", exp, ")"),
         x = "a_infl", y = "b_infl")
  ggsave(file.path(out_dir, paste0("influence_params_ab_", model_name, "_exp", exp, ".pdf")),
         plot = p, width = 6, height = 5)
}

plot_predictor_density <- function(df, pred_label, model_name, exp, out_dir) {
  p <- ggplot(df, aes(x = predictor)) +
    geom_density(fill = "steelblue", alpha = 0.4) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
    theme_minimal() +
    labs(title = paste0("Distribution of ", pred_label, " (", model_name, ", exp", exp, ")"),
         x = pred_label, y = "Density")
  ggsave(file.path(out_dir, paste0("predictor_density_", model_name, "_exp", exp, ".pdf")),
         plot = p, width = 7, height = 5)
}

## ---------- LOO ----------

compute_loo <- function(posterior_df, var_pattern = "^log_lik\\[") {
  vars <- names(posterior_df)[grepl(var_pattern, names(posterior_df))]
  mat  <- as.matrix(posterior_df[, vars])
  loo(mat, cores = parallel::detectCores())
}

plot_loo_comparison <- function(loo_list, title, filename) {
  comp <- loo_compare(loo_list)
  df_loo <- tibble(
    model        = rownames(comp),
    elpd_diff    = comp[, "elpd_diff"],
    se_elpd_diff = comp[, "se_diff"]
  )

  p <- ggplot(df_loo, aes(x = reorder(model, elpd_diff), y = elpd_diff, fill = model)) +
    geom_bar(stat = "identity", show.legend = FALSE) +
    geom_errorbar(aes(ymin = elpd_diff - se_elpd_diff, ymax = elpd_diff + se_elpd_diff),
                  width = 0.2, color = "black") +
    theme_minimal() +
    labs(title = title, x = "Model", y = "ELPD Difference")
  ggsave(filename, plot = p, width = 7, height = 5)
  comp
}

plot_loo_decomposed <- function(loo_total, loo_choice, loo_infl, title, filename) {
  fmt <- function(x, label) {
    tibble(
      model        = rownames(x),
      elpd_diff    = as.data.frame(x)[, "elpd_diff"],
      se_elpd_diff = as.data.frame(x)[, "se_diff"],
      component    = label
    )
  }
  df_all <- bind_rows(fmt(loo_total, "Total"), fmt(loo_choice, "Choice"), fmt(loo_infl, "Influence")) %>%
    mutate(component = factor(component, levels = c("Total", "Choice", "Influence")))

  p <- ggplot(df_all, aes(x = model, y = elpd_diff, fill = component)) +
    geom_bar(stat = "identity", position = position_dodge(width = 0.7), width = 0.6) +
    geom_errorbar(aes(ymin = elpd_diff - se_elpd_diff, ymax = elpd_diff + se_elpd_diff),
                  position = position_dodge(width = 0.7), width = 0.2, color = "black") +
    scale_fill_manual(values = c(Total = "grey30", Choice = "steelblue", Influence = "coral")) +
    theme_minimal() +
    labs(title = title, x = "Model", y = "ELPD Difference", fill = "")
  ggsave(filename, plot = p, width = 9, height = 5)
}
