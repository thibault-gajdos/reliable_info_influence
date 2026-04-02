###############################################################################
##  analysis_influence.r
##
##  Before running this script for the first time, run rename_lo_fit.r
##  to add named group-mean variables to the lo fit (avoids re-fitting).
##
##  This script:
##    1. Loads data (shared)
##    2. Loops over both models (lo & proba): diagnostics, PPCs, influence plots
##    3. Compares models via LOO (total, choice-only, influence-only)
###############################################################################

##rm(list = ls(all = TRUE))
source('~/thib/projects/tools/R_lib.r')
setwd('/Users/thibault/thib/projects/reliable_info/reliable_info_influence/')
source('analysis_influence_tools.r')

exp <- '7'

# ---------------------------------------------------
## MODEL SPECIFICATIONS
# ---------------------------------------------------
# var_prefix:  name of the Stan array containing the influence predictor
# pred_label:  human-readable label for the predictor in plot axes
# short_name:  short label used in plot titles and LOO legends
#
# lo / proba:       original models (all queried samples)
# lo_cc / proba_cc: choice-consistent variant (only queried samples
#                   matching the chosen color contribute to influence)

model_specs <- list(
  linear_choice_influence_lo = list(
    var_prefix = "sub_ev_out",
    pred_label = "sub_scaled",
    short_name = "LO"
  ),
  linear_choice_influence_lo_norm = list(
    var_prefix = "sub_ev_out",
    pred_label = "sub_scaled_norm",
    short_name = "LO_norm"
  ),
  linear_choice_influence_proba = list(
    var_prefix = "delta_p_out",
    pred_label = "delta_p",
    short_name = "Proba"
  ),
  linear_choice_influence_const_norm = list(
    var_prefix = "sub_ev_out",
    pred_label = "sub_scaled",
    short_name = "loo_const_norm"
  )
)




# ---------------------------------------------------
## LOAD AND PREPARE DATA (shared across models)
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

N       <- data_list$N
T_max   <- data_list$T_max
I_max   <- data_list$I_max
subjs   <- unique(data$ParticipantPrivateID)
t_subjs <- data_list$Tsubj



# ===================================================================
## PER-MODEL ANALYSIS LOOP
# ===================================================================
loo_total_list  <- list()
loo_choice_list <- list()
loo_infl_list   <- list()

for (model_name in names(model_specs)) {

out_dir <- paste0("./results/plots/exp", exp,'/',model_name)
sum_dir <- paste0("./results/summary/exp", exp,'/',model_name)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(sum_dir, recursive = TRUE, showWarnings = FALSE)

  cat(sprintf("\n\n========== %s ==========\n", model_name))
  spec       <- model_specs[[model_name]]
  short_name <- spec$short_name

  # --- Load posterior draws ---
  posterior_df <- load_posterior(model_name, exp)

  # --- Diagnostics ---
  plot_pairs(posterior_df, short_name, exp, out_dir)
  plot_traces(posterior_df, short_name, exp, out_dir)

  # --- Parameters ---
  gm           <- extract_group_means(posterior_df)
  params_indiv <- extract_indiv_params(posterior_df)

  cat("\n--- Group-level means ---\n")
  cat(sprintf("  w1..w5  = %.3f, %.3f, %.3f, %.3f, %.3f\n",
              gm$mu_w1, gm$mu_w2, gm$mu_w3, gm$mu_w4, gm$mu_w5))
  cat(sprintf("  alpha   = %.3f,  beta = %.3f\n", gm$mu_alpha, gm$mu_beta))
  cat(sprintf("  a_infl  = %.3f,  b_infl = %.3f,  sigma_infl = %.3f\n",
              gm$mu_a, gm$mu_b, gm$mu_sigma))

  # --- Summary tables ---
  save_summary(posterior_df, short_name, exp, sum_dir)

  # --- Choice plots ---
  plot_proba_transform(gm$mu_alpha, gm$mu_beta, params_indiv, short_name, exp, out_dir)
  plot_seq_weights(gm, params_indiv, short_name, exp, out_dir)
  plot_choice_calibration(posterior_df, data_list, t_subjs, N, short_name, exp, out_dir)
  plot_choice_calibration_per_subj(posterior_df, data_list, t_subjs, N, short_name, exp, out_dir)

  # --- Influence analysis ---
  df_infl <- build_infl_check_df(posterior_df, spec$var_prefix,
                                 data_list, data, subjs, t_subjs, N,
                                 params_indiv)
  cat(sprintf("  Influence check: %d rows\n", nrow(df_infl)))

  plot_infl_by_reliability(df_infl, short_name, exp, out_dir)
  plot_infl_by_reliability_per_subj(df_infl, short_name, exp, out_dir)
  plot_infl_demeaned_scatter(df_infl, spec$pred_label, short_name, exp, out_dir)
  plot_infl_within_cor(df_infl, spec$pred_label, short_name, exp, out_dir)
  plot_infl_slope_recovery(df_infl, spec$pred_label, short_name, exp, out_dir)
  plot_infl_demeaned_density(df_infl, posterior_df, short_name, exp, out_dir)
  plot_infl_by_proba_facet(df_infl, spec$pred_label, short_name, exp, out_dir)
  plot_infl_ab_scatter(gm, params_indiv, short_name, exp, out_dir)
  plot_predictor_density(df_infl, spec$pred_label, short_name, exp, out_dir)

  # --- LOO (store under short_name for readable legends) ---
  loo_total_list[[short_name]]  <- compute_loo(posterior_df, "^log_lik\\[")
  loo_choice_list[[short_name]] <- compute_loo(posterior_df, "^log_lik_choice\\[")
  loo_infl_list[[short_name]]   <- compute_loo(posterior_df, "^log_lik_influence\\[")

  cat(sprintf("  LOO elpd = %.1f\n",
              loo_total_list[[short_name]]$estimates["elpd_loo", "Estimate"]))
}

# ===================================================================
## LOO COMPARISON
# ===================================================================
cat("\n\n========== LOO COMPARISON ==========\n")

out_loo_dir <- paste0("./results/loo/exp", exp)
dir.create(out_loo_dir, recursive = TRUE, showWarnings = FALSE)

comp_total <- plot_loo_comparison(
  loo_total_list,
  paste0("LOO — total (exp", exp, ")"),
  file.path(out_loo_dir, paste0("loo_total_exp", exp, ".pdf"))
)
cat("\n--- Total ---\n"); print(comp_total)

comp_choice <- plot_loo_comparison(
  loo_choice_list,
  paste0("LOO — choice only (exp", exp, ")"),
  file.path(out_loo_dir, paste0("loo_choice_exp", exp, ".pdf"))
)
cat("\n--- Choice ---\n"); print(comp_choice)

comp_infl <- plot_loo_comparison(
  loo_infl_list,
  paste0("LOO — influence only (exp", exp, ")"),
  file.path(out_loo_dir, paste0("loo_influence_exp", exp, ".pdf"))
)
cat("\n--- Influence ---\n"); print(comp_infl)

plot_loo_decomposed(
  comp_total, comp_choice, comp_infl,
  paste0("LOO — decomposed (exp", exp, ")"),
  file.path(out_loo_dir, paste0("loo_decomposed_exp", exp, ".pdf"))
)

## Save comparison tables
write.csv(as.data.frame(comp_total),  file.path(out_loo_dir, paste0("loo_total_exp", exp, ".csv")))
write.csv(as.data.frame(comp_choice), file.path(out_loo_dir, paste0("loo_choice_exp", exp, ".csv")))
write.csv(as.data.frame(comp_infl),   file.path(out_loo_dir, paste0("loo_influence_exp", exp, ".csv")))

cat("\n=== Done. ===\n")
