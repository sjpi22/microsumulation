###########################  Benchmarking simulations  ################
#
#  Objective: Run benchmarking simulations
########################### <<<<<>>>>> ##############################################

rm(list = ls()) # Clean environment
options(scipen = 999) # View data without scientific notation

#### 1.Libraries and functions  ==================================================

###### 1.1 Load packages
library(tidyverse)
library(data.table)
library(ggplot2)

###### 1.2 Load functions

# Load functions
distr.sources <- list.files("R", 
                            pattern="*.R$", full.names=TRUE, 
                            ignore.case=TRUE, recursive = TRUE)
sapply(distr.sources, source, .GlobalEnv)

# Function for calculating prevalence with set arguments except method
prevalence_method <- function(method) {
  calc_prevalence(
    m_patients, 
    start_var = "time_P", 
    end_var = "time_C", 
    censor_var = "time_D", 
    method = method,
    v_ages = v_ages)
}

#### 2. General parameters ========================================================

###### 2.1 General parameters
n_cohort   <- 100000   # Number of individuals
p_cancer   <- 0.5     # Percentage at risk of preclinical cancer
seed       <- 123     # Random seed
conf_level <- 0.95    # Confidence level
n_sim      <- 1000    # Number of simulations
v_times    <- seq(0, 100, 0.5) # Time points for Weibull distribution
plt_size_text <- 18
plt_size_small <- 5/6*plt_size_text
path_weibull <- "benchmarking/distr_weibull.pdf" # Path to save Weibull distribution plot
path_benchmark <- "benchmarking/result.rds" # Path to save results

###### 2.2 Time-to-event parameters
# Individuals die independently of preclinical cancer, so effectively
# prevalence = CDF of preclinical Weibull distribution
l_params_cancer <- list(P  = list(distr = "weibull",
                                  params = list(shape = 4,
                                                scale = 60)), # Rate from birth to preclinical cancer onset
                        PC = list(distr = "unif",
                                  params = list(min = 150,
                                                max = 150)), # Rate from preclinical to clinical cancer
                        CD = list(distr = "unif",
                                  params = list(min = 150,
                                                max = 150)), # Rate from clinical cancer to death
                        Do = list(distr = "unif",
                                  params = list(min = 30,
                                                max = 100))) # Rate from birth to death from other causes # Rate from clinical cancer to death

# Parameters for individuals not at risk for cancer
l_params_norisk <- l_params_cancer
l_params_norisk$P <- l_params_cancer$PC

###### 2.3 Epidemiology calculation parameters
v_methods <- c("cs", "long", "rcs")
var_onset <- "time_P"
l_age_exp <- list(
  exp1 = c(30, 80),
  exp2 = seq(30, 80, 10) # Age ranges for prevalence
)


#### 3. Pre-processing  ===========================================

# Set seed for reproducibility
set.seed(seed)

# Create data to plot Weibull distribution
df_weibull <- data.frame(t = v_times,
                         cdf = query_distr("p", v_times, l_params_cancer[["P"]][["distr"]], l_params_cancer[["P"]][["params"]]))

# Plot Weibull distribution
plot_weibull <- ggplot(df_weibull, 
       aes(x = t, y = cdf)) +
  geom_line() +
  scale_fill_hue(h = c(180, 300)) +
  labs(x = "Age",
       y = "CDF") +
  theme_bw() + 
  theme(plot.title = element_blank(),
        axis.text.x = element_text(size = plt_size_small),
        axis.text.y = element_text(size = plt_size_small),
        axis.title.x = element_text(size = plt_size_small),
        axis.title.y = element_text(size = plt_size_small),
        legend.title = element_text(size = plt_size_small),
        legend.text = element_text(size = plt_size_small))

ggsave(path_weibull, 
       plot = plot_weibull, 
       width = 6, height = 4, dpi = 300)


#### 4. Analysis  ===========================================

# Initialize variables to store results
l_true_prevalence <- list()
final_time <- list()
final_prevalence <- list()

# Loop over age ranges to calculate true prevalence
for (i in names(l_age_exp)) {
  # Solve for true average prevalence in each age range
  v_ages <- l_age_exp[[i]]
  l_true_prevalence[[i]] <- mapply(function(a, b) {
    integrate(function(x){
      # Proportion with cancer weighted by proportion alive
      (1-query_distr("p", x, l_params_cancer[["Do"]][["distr"]], l_params_cancer[["Do"]][["params"]]))*
        query_distr("p", x, l_params_cancer[["P"]][["distr"]], l_params_cancer[["P"]][["params"]])}, a, b)$value /
      integrate(function(x){
        # Scale by weight alive
        1-query_distr("p", x, l_params_cancer[["Do"]][["distr"]], l_params_cancer[["Do"]][["params"]])
      }, a, b)$value
  }, 
  head(v_ages, -1), 
  v_ages[-1]) * p_cancer
}

# Initialize variables to store results
res_time <- list()
res_prevalence <- list()
res_summ <- list()
for (i in names(l_age_exp)) {
  # Time
  res_time[[i]] <- rep(0.0, length(v_methods))
  names(res_time[[i]]) <- v_methods
  
  # Prevalence summary stats
  res_prevalence[[i]] <- list()
  for (method in v_methods) {
    res_prevalence[[i]][[method]] <- data.table()
  }
  res_summ[[i]] <- data.table()
}

# Run simulations
for (i in 1:n_sim) {
  # Simulate cohort
  m_patients_cancer <- cancer_des(p_cancer*n_cohort, l_params_cancer) # Cancer cohort
  m_patients_norisk <- cancer_des((1-p_cancer)*n_cohort, l_params_norisk) # Non-cancer cohort
  m_patients <- rbind(m_patients_cancer, m_patients_norisk)
  
  # Loop over age ranges and methods
  for (i in names(l_age_exp)) {
    v_ages <- l_age_exp[[i]]
    for (method in v_methods) {
      # Calculate prevalence and time
      stime <- system.time({
        summ_prevalence <- prevalence_method(method)
      })
      
      # Save results and add time
      res_prevalence[[i]][[method]] <- rbind(res_prevalence[[i]][[method]], summ_prevalence)
      res_time[[i]][method] <- res_time[[i]][method] + stime["elapsed"]
    }
  }
}

# Summarize mean and SD by age and method
for (i in names(l_age_exp)) {
  for (method in v_methods) {
    dt_temp <- res_prevalence[[i]][[method]][, .(mean = mean(value),
                                      sd = sd(value)), 
                                  by = c("age_start", "age_end")][,
                                                                  `:=` (method = method,
                                                                        true = l_true_prevalence[[i]])]
    res_summ[[i]] <- rbind(res_summ[[i]], dt_temp)
  }
}

# Save results
saveRDS(list(time = res_time,
             outcomes = res_summ),
        file = path_benchmark)