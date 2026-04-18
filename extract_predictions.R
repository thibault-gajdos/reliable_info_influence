rm(list=ls(all=TRUE))  ## efface les données
##source('~/thib/projects/tools/R_lib.r')
##setwd('~/thib/projects/reliable_info')
##source('~/thib/projects/reliable_info/utils.r')

library(rstan)
library(tidyverse)
library(posterior)
library("writexl")
library(data.table)

# functions 
## Small function to compute the mode
# comp_mode <- function(x) {
#   ux <- unique(x)
#   ux[which.max(tabulate(match(x, ux)))]
# }
# 
# generate_pred <- function(fit){
#   #fit.name <- fit.name
#   pred_df = data.frame(subject = numeric(), trial = numeric(), pred = numeric())
#   #fit <- readRDS(fit.name)
#   pred <- extract(fit)$y_pred  ## (n_chains x Chains_length)x(n_subj)x(n_trials)
#   n_subj = dim(pred)[2]
#   for (i in c(1:n_subj)){
#     pred_temp <- extract(fit)$y_pred[,i,] 
#     pred_temp <-  as.data.frame(pred_temp) %>% 
#       summarize_all(comp_mode) %>%
#       t() %>%
#       as.data.frame() %>%
#       mutate(trial =row_number()) %>%
#       rename(pred = V1) %>%
#       mutate(subject = i) %>%
#       mutate(s = i) %>% ## subject index
#       mutate(trial = row_number())                
#     pred_df <- rbind(pred_df, pred_temp)
#   }
#   return(pred_df)
# }
# 

# load models
exp = '7';
model = c('fit_linear_choice_influence_proba')
load(paste0('./results/fits/exp',exp,'/',model,'_exp',exp,'.rdata') )
posterior_samples <- fit$draws()
posterior_df <- as_draws_df(posterior_samples)
selected_params <- posterior_df[, grepl("^pred_influence", colnames(posterior_df) ) ]
mean_vals <- sapply(selected_params, function(var) mean(var))
mean_df <- data.frame(
  parameter = names(mean_vals),
  mean = as.numeric(mean_vals)
)

write_xlsx(mean_df,  paste0('./results/predictions/exp',exp,'/influence_',model,'_exp',exp,'.xlsx') )

# 
# for (m in 1:4){
#   model = models[m]
#   fit <- readRDS(paste0('../../results/fits/exp',exp,'/',model,'_exp',exp,'.rds') )
#   
#   pred_df <- generate_pred(fit)
#   
#   write_xlsx(pred_df,  paste0('../../results/prediction/Exp',exp,'/ypred_',model,'_exp',exp,'.xlsx') )
# }
# 




