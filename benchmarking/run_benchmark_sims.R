###########################  Benchmarking computation time  ################
#
#  Objective: Run benchmarking simulations
########################### <<<<<>>>>> ##############################################

rm(list = ls()) # Clean environment
options(scipen = 999) # View data without scientific notation

#### 1.Libraries and functions  ==================================================

###### 1.1 Load packages
library(testthat)
library(tidyverse)
library(readxl)
library(data.table)
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
conf_level <- 0.95    # Confidence level
n_sim      <- 100     # Number of simulations

###### 2.2 Time-to-event parameters
l_params <- list(P  = list(distr = "weibull",
                           params = list(shape = 2,
                                         scale = 200)), # Rate from birth to preclinical cancer onset
                 PC = list(distr = "unif",
                           params = list(min = 150,
                                         max = 150)), # Rate from preclinical to clinical cancer
                 Do = list(distr = "unif",
                           params = list(min = 30,
                                         max = 100)), # Rate from birth to death from other causes
                 CD = list(distr = "unif",
                           params = list(min = 150,
                                         max = 150))) # Rate from clinical cancer to death

###### 2.3 Epidemiology calculation parameters
var_onset <- "time_P"
v_ages    <- seq(30, 80, 10) # Age ranges for prevalence


###### 4.5 Variation of prevalence
# Set seed for parallelization
set.seed(seed, kind = "L'Ecuyer-CMRG")

# If running locally, use all available cores except for reserved ones
registerDoParallel(cores = detectCores(logical = TRUE) - 2)

# Solve for true prevalence in each age range
true_prevalence <- df_cancer_cohort_wide[age_start < max(v_ages), 
                                         .(p_cases = sum(P),
                                           p_total = sum(H, P, C)),
                                         by = age_start]

# Run simulations
stime <- system.time({
  full_summ_prevalence <- foreach(
    i=1:n_sim, 
    .combine=rbind, 
    .inorder=FALSE, 
    .packages=c("data.table","tidyverse")) %dopar% {
      # Simulate cohort
      m_patients <- cancer_des(n_cohort, l_params)
      
      # Calculate prevalence (cross-sectional)
      summ_prevalence <- calc_prevalence(
        m_patients, 
        start_var = "time_P", 
        end_var = "time_C", 
        censor_var = "time_D", 
        method = "cs",
        v_ages = v_ages,
        output_uncertainty = T)
      
      # Calculate prevalence (longitudinal)
      summ_prevalence_long <- calc_prevalence(
        m_patients, 
        start_var = "time_P", 
        end_var = "time_C", 
        censor_var = "time_D", 
        method = "long",
        v_ages = v_ages,
        output_uncertainty = T)
      
      # Calculate prevalence (repeated cross-sectional)
      summ_prevalence_rcs <- calc_prevalence(
        m_patients, 
        start_var = "time_P", 
        end_var = "time_C", 
        censor_var = "time_D", 
        method = "rcs",
        v_ages = v_ages,
        output_uncertainty = T)
      
      # Merge longitudinal and cross-sectional
      summ_prevalence <- merge(summ_prevalence,
                               summ_prevalence_long, 
                               by = c("age_start", "age_end"),
                               suffixes = c("_cs", "_long"))
      
      # Merge with repeated cross-sectional
      summ_prevalence <- merge(summ_prevalence,
                               summ_prevalence_rcs, 
                               by = c("age_start", "age_end"),
                               suffixes = c("", "_rcs"))
      summ_prevalence
    }
})
print(stime)

# Calculate mean for each age group
mean_prevalence <- full_summ_prevalence[, .(mean_cs = mean(value_cs),
                                            mean_long = mean(value_long),
                                            mean_rcs = mean(value)), by = age_start]

# Calculate percentage bias
mean_prevalence[, `:=` (bias_cs = (mean_cs-v_true)/v_true,
                        bias_long = (mean_long-v_true)/v_true,
                        bias_rcs = (mean_rcs-v_true)/v_true)]

# Check whether CIs contains longitudinal and true prevalence
full_summ_prevalence[, `:=` (contained_true_cs = v_true >= ci_lb_cs & v_true <= ci_ub_cs,
                             contained_true_long = v_true >= ci_lb_long & v_true <= ci_ub_long,
                             contained_true_rcs = v_true >= ci_lb & v_true <= ci_ub,
                             consistent_cs_long = value_long >= ci_lb_cs & value_long <= ci_ub_cs,
                             consistent_cs_rcs = value >= ci_lb_cs & value <= ci_ub_cs)]

# Percentage of values within CIs
pct_contained <- full_summ_prevalence[, .(pct_cs = mean(contained_true_cs),
                                          pct_long = mean(contained_true_long),
                                          pct_rcs = mean(contained_true_rcs),
                                          pct_consistent_long = mean(consistent_cs_long),
                                          pct_consistent_rcs = mean(consistent_cs_rcs)), by = age_start]

# Perform consistency unit tests
test_that("Methods of calculating prevalence match with confidence intervals", {
  expect_equal(abs(pct_contained$pct_cs - conf_level) < 0.03, rep(T, length(v_ages)-1))
  # expect_equal(abs(pct_contained$pct_long - conf_level) < 0.03, rep(T, length(v_ages)-1))
  # expect_equal(abs(pct_contained$pct_rcs - conf_level) < 0.03, rep(T, length(v_ages)-1))
})

test_that("Methods of calculating prevalence produce similar results", {
  expect_equal(pct_contained$pct_consistent_long > conf_level, rep(T, length(v_ages)-1))
  expect_equal(pct_contained$pct_consistent_rcs > conf_level, rep(T, length(v_ages)-1))
})
