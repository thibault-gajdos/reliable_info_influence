rm(list=ls(all=TRUE))  ## efface les données
source('~/thib/projects/tools/R_lib.r')
setwd('/Users/thibault/thib/projects/reliable_info/reliable_info_identif/')

##################################################
##        PREPARE THE DATA
##################################################

exp <- 7
# specify which experiment to run
## target <- 'proba'
load(paste0("./data/data_list_exp", as.character(exp), ".rdata"))
# data_list <- data_list_choice
# data_list <- data_list_proba


data_list$grainsize = 5 ## specify grainsize for within chain parallelization

#####################################################
##  Fit THE MODELS
####################################################

model_name <- 'linear_choice_influence'
## Compile the model
model <- cmdstan_model(
    stan_file = paste0('./stan/',model_name,'.stan'),
    force_recompile = TRUE, ## necessary if you change the mode
    cpp_options = list(stan_opencl = FALSE, stan_threads = TRUE), ## within chain parallel
    stanc_options = list("O1"), ## fastest sampling
    compile_model_methods = TRUE ## necessary for loo moment matching
)

##MAP
## fit_map <- model$optimize(
##         data = data_list,
##         ##init = fit_map$draws(),
##         init = 0,
##         algorithm = "lbfgs",
##         threads = 5
## )
## mle <- fit_map$mle()
## mle["mu_theta_choice"]
## mle["mu_kappa"]
## mle["mu_Psi"]
## mle["mu_l_anchor"]
## mle["mu_Delta_plus"]
## mle["mu_Delta_minus"]
## mle["mu_lambda_recency"]
## log_lik_vec <- mle[grep("^log_lik\\[", names(mle))]
## total_loglik <- sum(log_lik_vec)
## total_loglik

## mle["mu_alpha"]
## mle["mu_beta"]
## mle["mu_lambda"]
## mle["mu_theta"]


## Sampling
dir.create(paste0("./results/fits/exp", as.character(exp), "/CSV_", model_name, "BLO/"), showWarnings =   FALSE)

fit <- model$sample(
  data = data_list,
  seed = 4321,
  ##init = list(inits_chain,inits_chain, inits_chain,inits_chain),
  ##init =0,
  ##init = fit_map,
  chains = 4,
  parallel_chains = 4,
  threads_per_chain = 5,
  iter_warmup = 3000,
  iter_sampling = 1000,
  max_treedepth = 12,
  adapt_delta = .95,
  save_warmup = FALSE
)


## Compute LOO
loo <- fit$loo(cores = 1, moment_match = TRUE)

## Save results
dir.create(paste0('./results/fits/exp', as.character(exp),'/'), showWarnings = FALSE)
save(fit, file = paste0("./results/fits/exp", as.character(exp), "/fit_", model_name, "_exp", as.character(exp), ".rdata"))
fit$save_object(paste0("./results/fits/exp",as.character(exp),"/fit_",model_name,"_exp",as.character(exp),".rds"))

dir.create(paste0('./results/loo/exp', as.character(exp),'/'), showWarnings = FALSE)
save(loo, file = paste0("./results/loo/exp",as.character(exp),"/loo_",model_name,"_exp",as.character(exp),".rdata"))
## to load:
## fit <- cmdstanr::read_cmdstanr("./results/fits/exp8/fit_BLO_choice.rds")

