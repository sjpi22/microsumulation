###########################  Unit test: Risk   ################
#
#  Objective: Run unit tests for calculating age-conditional and lifetime
########################### <<<<<>>>>> ##############################################

rm(list = ls()) # Clean environment
options(scipen = 999) # View data without scientific notation

#### 1.Libraries and functions  ==================================================

###### 1.1 Load packages
library(testthat)
library(tidyverse)
library(readxl)
library(data.table)
require(deSolve)
library(foreach)
library(doParallel)

###### 1.1 Load functions

# Load functions
distr.sources <- list.files("R", 
                            pattern="*.R$", full.names=TRUE, 
                            ignore.case=TRUE, recursive = TRUE)
sapply(distr.sources, source, .GlobalEnv)


#### 2. General parameters ========================================================

###### 2.1 General parameters
n_cohort   <- 50000   # Number of individuals
max_age    <- 100     # Maximum age
seed       <- 123     # Random seed
n_sim      <- 500     # Number of simulations
conf_level <- 0.95

###### 2.2 Time-to-event parameters
l_params <- list(r_P  = 1/200, # Rate from birth to preclinical cancer onset
                 r_PC = 1/10, # Rate from preclinical to clinical cancer
                 r_Do = 1/80, # Rate from birth to death from other causes
                 r_CD = 1/10) # Rate from clinical cancer to death

###### 2.3 Epidemiology calculation parameters
var_onset <- "time_P"
v_ages    <- c(30, 40, 50) # Age ranges for incidence


#### 3. Pre-processing  ===========================================

# Set seed for reproducibility
set.seed(seed)


#### 4. Unit tests  ===========================================

###### 4.1 Lifetime risk
# Set seed for parallelization
set.seed(seed, kind = "L'Ecuyer-CMRG")

# If running locally, use all available cores except for reserved ones
registerDoParallel(cores = detectCores(logical = TRUE) - 2)

# Vector with initial states
v_state_init <- c(H  = 1, 
                  P  = 0, 
                  C  = 0,
                  DO = 0, 
                  DC = 0,
                  CInc = 0)

# Solves the system of ODEs an returns the proportion or number of the 
# population in each of the states or compartments at the user-specified times
# in a data.frame in wide format
df_cancer_cohort_wide <- as.data.table(lsoda(y     = v_state_init, 
                                             times = 0:max_age, 
                                             func  = cancer_cohort_ode, 
                                             parms = l_params))

# Calculate true risk for scenarios
v_true <- c(
  df_cancer_cohort_wide[time == max_age, CInc],
  df_cancer_cohort_wide[time == max(v_ages), CInc],
  (df_cancer_cohort_wide[time == max(v_ages), CInc] - 
     df_cancer_cohort_wide[time == min(v_ages), CInc])/
    df_cancer_cohort_wide[time == min(v_ages), H + P]
)

# Run simulations
stime <- system.time({
  full_summ_incidence <- foreach(
    i=1:n_sim, 
    .combine=rbind, 
    .inorder=FALSE, 
    .packages=c("data.table","tidyverse")) %dopar% {
      # Simulate cohort
      m_patients <- cancer_des(n_cohort, l_params)
      
      summ_risk <- rbind(
        # Calculate lifetime risk
        calc_risk(
          m_patients, 
          start_var = "time_C", 
          censor_var = "time_D", 
          min_age = 0,
          max_age = max_age,
          output_uncertainty = T) %>%
          mutate(label = "lifetime"),
        # Calculate upper bounded risk
        calc_risk(
          m_patients, 
          start_var = "time_C", 
          censor_var = "time_D",
          max_age = max(v_ages),
          output_uncertainty = T) %>%
          mutate(label = "ub"),
        # Calculate risk in an age range
        calc_risk(
          m_patients, 
          start_var = "time_C", 
          censor_var = "time_D", 
          min_age = min(v_ages),
          max_age = max(v_ages),
          output_uncertainty = T) %>%
          mutate(label = "range")
      )
      summ_risk
    }
})
print(stime)

# Calculate mean for each age group
mean_risk <- full_summ_incidence[, .(mean_value = mean(value)), by = label]

# Bias
mean_risk[, `:=` (bias = (mean_value - v_true)/v_true)]

# Check whether CIs contains longitudinal and true incidence
full_summ_incidence[, `:=` (contained_true = v_true >= ci_lb & v_true <= ci_ub)]

# Percentage of values within CIs
pct_contained <- full_summ_incidence[, .(pct = mean(contained_true)), by = label]

# Calculate SD for each age group
sd_incidence <- full_summ_incidence[, .(mcse = sd(value),
                                        mean_se = mean(se)), by = label]

# Calculate ratios of SDs
sd_incidence[, `:=` (ratio_long = mcse / mean_se)]

# Perform consistency unit tests
test_that("Methods of calculating incidence produce similar results", {
  expect_equal(abs(pct_contained$pct - conf_level) < 0.03, rep(T, 3))
})
