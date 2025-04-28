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

# Function for calculating incidence with set arguments except method
incidence_method <- function() {
  calc_incidence(
    m_patients, 
    time_var = "time_P", 
    censor_var = "time_D", 
    rate_unit = rate_unit,
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
path_output <- "benchmarking/benchmark_output.rds" # Path to save output

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
rate_unit <- 100000 # Incidence rate unit


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

# Loop over age ranges to calculate true prevalence and incidence
l_true_prevalence <- list()
l_true_incidence <- list()
for (i in names(l_age_exp)) {
  # Extract age range
  v_ages <- l_age_exp[[i]]
  
  # Solve for true average prevalence in each age range
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
  
  # Solve for true incidence in each age group
  l_true_incidence[[i]] <- mapply(function(a, b) {
    integrate(function(x){
      # Total proportion of entire population developing condition in time range:
      # Integrate over proportion alive times probability density of developing condition
      (1-query_distr("p", x, l_params_cancer[["Do"]][["distr"]], l_params_cancer[["Do"]][["params"]])) *
        query_distr("d", x, l_params_cancer[["P"]][["distr"]], l_params_cancer[["P"]][["params"]])
    }, a, b)$value /
      # Divide by living person-years in time period
      integrate(function(x){
        1-query_distr("p", x, l_params_cancer[["Do"]][["distr"]], l_params_cancer[["Do"]][["params"]])
      }, a, b)$value
  }, 
  head(v_ages, -1), 
  v_ages[-1]) * p_cancer * rate_unit
}

# Initialize variables to store results
full_time <- list()
full_stats <- list()
for (i in names(l_age_exp)) {
  # Prevalence summary stats
  full_stats[[i]] <- list()
  for (method in c(v_methods, "incidence")) {
    full_stats[[i]][[method]] <- data.table()
    full_time[[i]][[method]] <- c()
  }
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
      full_stats[[i]][[method]] <- rbind(full_stats[[i]][[method]], summ_prevalence)
      full_time[[i]][[method]] <- c(full_time[[i]][[method]], stime["elapsed"])
    }
    
    # Calculate incidence and time
    stime <- system.time({
      summ_incidence <- incidence_method()
    })
    
    # Save results and add time
    full_stats[[i]][["incidence"]] <- rbind(full_stats[[i]][["incidence"]], summ_incidence)
    full_time[[i]][["incidence"]] <- cbind(full_time[[i]][["incidence"]], stime["elapsed"])
  }
}

# Summarize mean and SD by age and method
res_summ <- list()
res_time <- list()
for (i in names(l_age_exp)) {
  res_summ[[i]] <- data.table()
  res_time[[i]] <- list()
  for (method in c(v_methods, "incidence")) {
    if (method == "incidence") {
      v_true <- l_true_incidence[[i]]
    } else {
      v_true <- l_true_prevalence[[i]]
    }
    
    # Summarize mean value
    dt_stats <- full_stats[[i]][[method]][, .(
      mean = mean(value),
      sd = sd(value)), 
      by = c("age_start", "age_end")][,
                                      `:=` (
                                        method = method,
                                        true = v_true)]
    res_summ[[i]] <- rbind(res_summ[[i]], dt_stats)
    
    # Summarize time
    res_time[[i]][[method]] <- c(mean = mean(full_time[[i]][[method]]), 
                                 sd = sd(full_time[[i]][[method]]))
  }
}

# Save results
saveRDS(list(time = res_time,
             outcomes = res_summ),
        file = path_output)
