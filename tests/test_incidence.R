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

###### 1.2 Load functions

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
conf_level <- 0.95    # Confidence level
n_sim      <- 500     # Number of simulations

###### 2.2 Time-to-event parameters
l_params_ode <- list(r_P  = 1/200, # Rate from birth to preclinical cancer onset
                     r_PC = 1/10, # Rate from preclinical to clinical cancer
                     r_CD = 1/10, # Rate from clinical cancer to death
                     r_Do = 1/80) # Rate from birth to death from other causes

###### 2.3 Epidemiology calculation parameters
var_onset <- "time_P"
v_ages    <- c(30, 40, 50) # Age ranges for incidence
rate_unit <- 1000 # Unit for incidence rate


#### 3. Pre-processing  ===========================================

# Set seed for reproducibility
set.seed(seed)


#### 4. Unit tests  ===========================================

###### 4.1 True incidence
# Set constant incidence
incidence_exp <- 0.2

# Initialize matrix of patient trajectories
m_patients <- data.table(pt_id = 1:6)
setkey(m_patients, pt_id)

# Set time to cancer and death
m_patients[, `:=` (time_C = c(NA, 18, 120, 20, 42, 80),
                   time_D = c(80, 20, 40, 45, 45, 85))]

# Calculate censor time
m_patients[, `:=` (time_censor = pmin(time_C, time_D, na.rm = T))]

# Test accuracy of longitudinal incidence
test_that("Test accuracy of longitudinal incidence", {
  # Calculate incidence out of entire population
  summ_incidence <- calc_incidence(
    m_patients, 
    time_var = "time_C", 
    censor_var = "time_D",
    v_ages = v_ages, 
    rate_unit = rate_unit)
  
  # Check expected value of incidence
  expect_equal(summ_incidence$value, c(0, 1/30)*rate_unit)
  
  # Check expected person-years with cancer
  expect_equal(summ_incidence$n_events, c(0, 1))
  
  # Check expected total person-years
  expect_equal(summ_incidence$person_years_total, c(50, 30))
  
  # Calculate incidence out of cancer-free population
  summ_incidence <- calc_incidence(
    m_patients, 
    time_var = "time_C", 
    censor_var = "time_censor",
    v_ages = v_ages, 
    rate_unit = rate_unit)
  
  # Check expected value of incidence
  expect_equal(summ_incidence$value, c(0, 1/22)*rate_unit)
  
  # Check expected person-years with cancer
  expect_equal(summ_incidence$n_events, c(0, 1))
  
  # Check expected total person-years
  expect_equal(summ_incidence$person_years_total, c(40, 22))
})


###### 4.2 Variation of incidence
# Vector with initial states
v_state_init <- c(H  = 1, 
                  P  = 0, 
                  C  = 0,
                  DO = 0, 
                  DC = 0,
                  CInc = 0)

# Generate DES parameters
l_params_des <- list()
for (event in names(l_params_ode)) {
  l_params_des[[substr(event, 3, nchar(event))]] <- list(
    distr = "exp",
    params = list(rate = l_params_ode[[event]]))
}

# Set seed for parallelization
set.seed(seed, kind = "L'Ecuyer-CMRG")

# If running locally, use all available cores except for reserved ones
registerDoParallel(cores = detectCores(logical = TRUE) - 2)

# Solves the system of ODEs an returns the proportion or number of the 
# population in each of the states or compartments at the user-specified times
# in a data.frame in wide format
df_cancer_cohort_wide <- as.data.table(lsoda(y     = v_state_init, 
                                             times = 0:max_age, 
                                             func  = cancer_cohort_ode, 
                                             parms = l_params_ode))

# Lag age and calculate incidence with diff
df_cancer_cohort_wide[, `:=` (time = time + 1,
                              p_alive = H + P + C,
                              p_cf = H + P,
                              dCInc = c(diff(CInc), NA))]

# Map age range groups and calculate total alive per age
df_cancer_cohort_wide[, `:=` (age_idx = findInterval(time, v_ages),
                              p_alive_avg = (p_alive + lead(p_alive))/2,
                              p_cf_avg = (p_cf + lead(p_cf))/2)]
df_cancer_cohort_wide[, `:=` (age_start = v_ages[age_idx]), by = age_idx]

# Calculate true incidence in age ranges
true_incidence <- df_cancer_cohort_wide[age_start < max(v_ages), 
                                        .(true = sum(dCInc)/sum(p_alive_avg)*rate_unit,
                                          true_cf = sum(dCInc)/sum(p_cf_avg)*rate_unit), 
                                        by = age_start]
v_true <- true_incidence$true
v_true_cf <- true_incidence$true_cf

# Run simulations
stime <- system.time({
  full_summ_incidence <- foreach(
    i=1:n_sim, 
    .combine=rbind, 
    .inorder=FALSE, 
    .packages=c("data.table","tidyverse")) %dopar% {
      # Simulate cohort
      m_patients <- cancer_des(n_cohort, l_params_des)
      
      # Calculate incidence (longitudinal, full)
      summ_incidence_long <- calc_incidence(
        m_patients, 
        time_var = "time_C", 
        censor_var = "time_D", 
        v_ages = v_ages,
        output_uncertainty = T,
        rate_unit = rate_unit)
      
      # Calculate incidence (longitudinal, cancer-free)
      summ_incidence_long_cf <- calc_incidence(
        m_patients, 
        time_var = "time_C", 
        censor_var = "time_C", 
        v_ages = v_ages,
        output_uncertainty = T,
        rate_unit = rate_unit)
      
      # Merge data
      summ_incidence_long <- merge(summ_incidence_long,
                                   summ_incidence_long_cf,
                                   by = c("age_range", "age_start", "age_end"),
                                   suffixes = c("", "_cf"))
      
      summ_incidence_long
    }
})
print(stime)

# Calculate mean for each age group
mean_incidence <- full_summ_incidence[, .(mean_long = mean(value),
                                          mean_long_cf = mean(value_cf)), by = age_start]

# Bias
mean_incidence[, `:=` (bias_long = (mean_long-v_true)/v_true,
                       bias_long_cf = (mean_long_cf-v_true_cf)/v_true_cf)]

# Check whether CIs contains longitudinal and true incidence
full_summ_incidence[, `:=` (contained_true_long = v_true >= ci_lb & v_true <= ci_ub,
                            contained_true_long_cf = v_true_cf >= ci_lb_cf & v_true_cf <= ci_ub_cf)]

# Percentage of values within CIs
pct_contained <- full_summ_incidence[, .(pct_long = mean(contained_true_long),
                                         pct_long_cf = mean(contained_true_long_cf)), by = age_start]

# Calculate SD for each age group
sd_incidence <- full_summ_incidence[, .(mcse_long = sd(value),
                                        mean_se_long = mean(se)), by = age_start]

# Calculate ratios of SDs
sd_incidence[, `:=` (ratio_long = mcse_long / mean_se_long)]

# Perform consistency unit tests
test_that("Methods of calculating incidence produce similar results", {
  expect_equal(abs(pct_contained$pct_long - conf_level) < 0.03, rep(T, length(v_ages)-1))
  expect_equal(abs(pct_contained$pct_long_cf - conf_level) < 0.03, rep(T, length(v_ages)-1))
})
