rm(list=ls(all=TRUE))  ## efface les données
##source('~/thib/projects/tools/R_lib.r')
##setwd('~/thib/projects/reliable_info')
##source('~/thib/projects/reliable_info/utils.r')

library(rstan)
library(tidyverse)
library(posterior)
library("writexl")
library(data.table)


# load models
exp = '7';
model = c('fit_linear_choice_influence_proba')
load(paste0('./results/fits/exp',exp,'/',model,'_exp',exp,'.rdata') )
posterior_samples <- fit$draws()
posterior_df <- as_draws_df(posterior_samples)

# extracting predictions of the influence reports
selected_params <- posterior_df[, grepl("^pred_influence", colnames(posterior_df) ) ]
mean_vals <- sapply(selected_params, function(var) mean(var))
mean_df <- data.frame(
  parameter = names(mean_vals),
  mean = as.numeric(mean_vals)
)
write_xlsx(mean_df,  paste0('./results/predictions/exp',exp,'/influence_',model,'_exp',exp,'.xlsx') )


# for (m in 1:4){
#   model = models[m]
#   fit <- readRDS(paste0('../../results/fits/exp',exp,'/',model,'_exp',exp,'.rds') )
#   
#   pred_df <- generate_pred(fit)
#   
#   write_xlsx(pred_df,  paste0('../../results/prediction/Exp',exp,'/ypred_',model,'_exp',exp,'.xlsx') )
# }
# 

# extracting group-level parameters
extract_group_summary <- function(posterior_df, model_name) {
  
  params <- c(
    "mu_alpha", "mu_beta", "mu_w1", "mu_w2", "mu_w3", "mu_w4", "mu_w5",
    "mu_a_infl", "mu_b_infl", "mu_sigma_infl"
  )
  
  get_stats <- function(name) {
    if (name %in% names(posterior_df)) {
      x <- posterior_df[[name]]
      c(
        mean   = mean(x, na.rm = TRUE),
        `2.5%`  = unname(quantile(x, 0.025, na.rm = TRUE)),
        `97.5%` = unname(quantile(x, 0.975, na.rm = TRUE))
      )
    } else {
      c(mean = NA_real_, `2.5%` = NA_real_, `97.5%` = NA_real_)
    }
  }
  
  mat <- t(sapply(params, get_stats))
  
  result <- data.frame(
    row       = seq_along(params),
    parameter = params,
    mean      = mat[, "mean"],
    `2.5%`    = mat[, "2.5%"],
    `97.5%`   = mat[, "97.5%"],
    model     = model_name,
    row.names = NULL,
    check.names = FALSE   # <-- this is the key fix
  )
  
  return(result)
}
mname = 'influence_proba'; 
group           <- extract_group_summary(posterior_df, mname)
write.csv(group, file = paste0('./results/summary/exp',exp,'/param_group_influence_',model,'_exp',exp,'.csv'))


# extracting individual-level parameters

extract_indiv_summary <- function(posterior_df, model_name) {
  
  params_vars <- names(posterior_df)[grepl("^params\\[[0-9]+,[0-9]+\\]$", names(posterior_df))]
  
  get_stats <- function(x) {
    c(
      mean   = mean(x, na.rm = TRUE),
      `2.5%`  = unname(quantile(x, 0.025, na.rm = TRUE)),
      `97.5%` = unname(quantile(x, 0.975, na.rm = TRUE))
    )
  }
  
  stats_mat <- t(sapply(params_vars, function(v) get_stats(posterior_df[[v]])))
  
  idx <- stringr::str_match(params_vars, "^params\\[([0-9]+),([0-9]+)\\]$")
  
  df <- data.frame(
    subj = as.integer(idx[, 2]),
    k    = as.integer(idx[, 3]),
    stats_mat,
    row.names = NULL,
    check.names = FALSE
  )
  
  # original mapping from k → parameter name
  param_names <- c(
    "w1", "w2", "w3", "w4", "w5",
    "alpha", "beta", "a_inf", "b_inf", "sigma_inf"
  )
  
  df$parameter <- param_names[df$k]
  df$model <- model_name
  
  # your desired order
  desired_order <- c(
    "alpha", "beta",
    "w1", "w2", "w3", "w4", "w5",
    "a_inf", "b_inf", "sigma_inf"
  )
  
  df$parameter <- factor(df$parameter, levels = desired_order)
  
  # sort by parameter, then subject
  df <- df[order(df$parameter, df$subj), ]
  
  return(df)
}

indv         <- extract_indiv_summary(posterior_df, mname)
write.csv(indv, file = paste0('./results/summary/exp',exp,'/param_individual_influence_',model,'_exp',exp,'.csv'))

## extract full distributions 

extract_group_draws <- function(posterior_df) {
  
  params <- c(
    "mu_w1", "mu_w2", "mu_w3", "mu_w4", "mu_w5",
    "mu_alpha", "mu_beta", "mu_a_infl", "mu_b_infl", "mu_sigma_infl"
  )
  
  # check all params exist
  missing <- setdiff(params, names(posterior_df))
  if (length(missing) > 0) {
    stop(paste("Missing parameters:", paste(missing, collapse = ", ")))
  }
  
  # extract as matrix (M samples × N parameters)
  mat <- as.matrix(posterior_df[, params])
  
  return(mat)
}

draws         <- extract_group_draws(posterior_df)
write.csv(draws, file = paste0('./results/summary/exp',exp,'/full_distribution_influence_',model,'_exp',exp,'.csv'))

