rm(list=ls(all=TRUE))  ## efface les données
source('~/thib/projects/tools/R_lib.r')
setwd('/Users/thibault/thib/projects/reliable_info/reliable_info_influence/')


##################################################
##        PREPARE THE DATA       
##################################################
exp <- 6
load(file = paste0("./data/data_reliability_exp", as.character(exp), ".rdata"))

##  if Response  ResponseButtonOrder= 1:  blue->1, red->0
##  if Response  ResponseButtonOrder= 0:  blue->0, red->1
##  we recode: blue ->1, red ->2
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
    mutate(influence = influence/100 - .5)

## influence_i = 1 if proba_i = influence_proba, and 0 otherwise
data <- data %>%
  mutate(across(proba_1:proba_6, ~ as.integer(.x == influence_proba), .names = "{.col}_inf")) %>%
  rename_with(~ str_replace(.x, "proba_(\\d+)_inf", "influence_\\1"))


N = length(unique(data$ParticipantPrivateID))
T_max = max(data$TrialNumber)
I_max <- max(data$sample_number) ## max number of samples/trial
## compute trials by subject
d <- data %>%
  group_by(ParticipantPrivateID) %>%
  summarise(t_subjs = n())
t_subjs <- d$t_subjs
subjs <- unique(data$ParticipantPrivateID)


## Initialize data arrays
choice <- array(-1, c(N, T_max))
influence_sample <- array(-1, c(N, T_max, I_max))
influence <- array(-1, c(N, T_max))
color <- array( -1, c(N, T_max, I_max))
proba <- array(-1, c(N, T_max, I_max))
sample <- array(-1, c(N, T_max))
## fill the  arrays
for (n in 1:N) { ## loop through subjects
  t <- t_subjs[n] ## number of trials for subj i
  data_subj <- data %>% filter(ParticipantPrivateID == subjs[n])
  choice[n, 1:t] <- data_subj$choice
  influence[n, 1:t] <- data_subj$influence
  for (k in 1:t) { ## loop through trials
    data_subj_t <- data_subj[k,]
    sample[n,k] <- data_subj_t$sample_number
    for (i in 1:data_subj_t$sample_number) {
      color_var <- paste0("color_", i)
      proba_var <- paste0("proba_", i)
      influence_var <- paste0("influence_", i)
      color[n, k, i] <- data_subj[[color_var]][k]
      proba[n, k, i] <- data_subj[[proba_var]][k] / 100
      influence_sample[n, k, i] <- data_subj[[influence_var]][k]
    }
  }
}


data_list <- list(
  N = N,
  T_max = T_max,
  I_max = I_max,
  Tsubj = t_subjs,
  color = color,
  proba = proba,
  choice = choice, ## 
  influence = influence, 
  influence_sample = influence_sample,
  sample = sample
)


save(data_list, file = paste0('./data/data_list_exp', as.character(exp), '.rdata'))

