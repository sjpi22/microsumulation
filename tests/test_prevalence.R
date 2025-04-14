###########################  Unit test: Prevalence   ################
#
#  Objective: Run unit tests for calculating prevalence
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

# Function to simulate data
simulate_data <- function(n, l_params) {
  # Initialize matrix of patient trajectories
  m_patients <- data.table(pt_id = 1:n)
  setkey(m_patients, pt_id)
  
  # Simulate time to death from other causes, preclinical cancer, clinical cancer, and death from cancer
  with(l_params, {
    m_patients[, `:=` (time_Do = runif(.N, min = params_Do["min"], max = params_Do["max"]),
                       time_P  = runif(.N, min = params_P["min"], max = params_P["max"]),
                       time_PC = runif(.N, min = params_PC["min"], max = params_PC["max"]),
                       time_CD = runif(.N, min = params_CD["min"], max = params_CD["max"]))]
    
    # Calculate time from birth to clinical cancer and death from cancer
    m_patients[, `:=` (time_C  = time_P + time_PC,
                       time_Dc = time_P + time_PC + time_CD)]
    
    # Calculate time to death and cause of death
    m_patients[, `:=` (time_D = pmin(time_Do, time_Dc),
                       Dc = (time_Dc < time_Do))]
  })
  
  return(m_patients)
}


#### 2. General parameters ========================================================

###### 2.1 General parameters
n_cohort   <- 50000   # Number of individuals
max_age    <- 100     # Maximum age
seed       <- 123     # Random seed
conf_level <- 0.95    # Confidence level
n_sim      <- 500     # Number of simulations

###### 2.2 Time-to-event parameters
params_Do <- c(min = 30, max = 100)  # Time from birth to death from other causes
params_P  <- c(min = 0, max = 500)  # Time from birth to preclinical cancer onset
params_PC <- c(min = 0, max = 10)  # Time from preclinical to clinical cancer
params_CD <- c(min = 0, max = 10)  # Time from clinical cancer to death
l_params <- list(params_Do = params_Do,
                 params_P = params_P,
                 params_PC = params_PC,
                 params_CD = params_CD)

###### 2.3 Epidemiology calculation parameters
var_onset <- "time_P"
v_ages    <- seq(30, 80, 10) # Age ranges for prevalence


#### 3. Pre-processing  ===========================================

# Set seed for reproducibility
set.seed(seed)


#### 4. Prevalence  ===========================================

###### 4.1 True prevalence - condition ends before death
# Set constant prevalence
prevalence_exp <- 0.2

# Initialize matrix of patient trajectories
m_patients <- data.table(pt_id = 1:n_cohort)
setkey(m_patients, pt_id)

# Simulate time to death and preclinical cancer at fixed intervals
m_patients[, `:=` (time_D = max_age,
                   time_P = rep(seq(0, max_age - max_age*prevalence_exp, by = max_age*prevalence_exp), 
                                each = n_cohort*prevalence_exp))]

# Calculate time from to clinical cancer as fixed interval from preclinical cancer
m_patients[, `:=` (time_C = time_P + max_age/5)]

# Test accuracy of cross-sectional prevalence
test_that("Test accuracy of cross-sectional prevalence", {
  summ_prevalence <- calc_prevalence(
    m_patients, 
    start_var = "time_P", 
    end_var = "time_C", 
    censor_var = "time_D", 
    method = "cs",
    v_ages = v_ages[2:3])
  
  expect_equal(summ_prevalence$value, prevalence_exp)
  expect_equal(summ_prevalence$n_cases, n_cohort*prevalence_exp)
})

# Test accuracy of longitudinal prevalence
test_that("Test accuracy of longitudinal prevalence", {
  summ_prevalence <- calc_prevalence(
    m_patients, 
    start_var = "time_P", 
    end_var = "time_C", 
    censor_var = "time_D", 
    method = "long",
    v_ages = v_ages)
  
  expect_equal(summ_prevalence$value, rep(prevalence_exp, length(v_ages)-1))
  expect_equal(summ_prevalence$person_years_cases, rep(n_cohort*prevalence_exp*diff(v_ages[1:2]), length(v_ages)-1))
  expect_equal(summ_prevalence$person_years_total, rep(n_cohort*diff(v_ages[1:2]), length(v_ages)-1))
})

###### 4.2 True prevalence - condition ends after death
# Set constant prevalence
prevalence_exp <- 0.5

# Initialize matrix of patient trajectories
m_patients <- data.table(pt_id = 1:n_cohort)
setkey(m_patients, pt_id)

# Simulate time to death and preclinical cancer at fixed intervals
m_patients[, `:=` (time_D = max_age,
                   time_P = rep(c(min(tail(v_ages, 2)), max_age + 1), each = n_cohort*prevalence_exp))]

# Calculate time from to clinical cancer as fixed interval from preclinical cancer
m_patients[, `:=` (time_C = time_P + 20)]

# Test accuracy of cross-sectional prevalence
test_that("Test accuracy of cross-sectional prevalence", {
  # Prevalence for early age should be 0
  summ_prevalence <- calc_prevalence(
    m_patients, 
    start_var = "time_P", 
    end_var = "time_C", 
    censor_var = "time_D", 
    method = "cs",
    v_ages = v_ages[1:2])
  
  expect_equal(summ_prevalence$value, 0)
  
  # Prevalence for last age category should be ~1/4
  summ_prevalence <- calc_prevalence(
    m_patients, 
    start_var = "time_P", 
    end_var = "time_C", 
    censor_var = "time_D", 
    method = "cs",
    v_ages = tail(v_ages, 2))
  
  expect_equal(summ_prevalence$value, prevalence_exp)
  expect_equal(summ_prevalence$n_cases, n_cohort*prevalence_exp)
})

# Test accuracy of longitudinal prevalence
test_that("Test accuracy of longitudinal prevalence", {
  summ_prevalence <- calc_prevalence(
    m_patients, 
    start_var = "time_P", 
    end_var = "time_C", 
    censor_var = "time_D", 
    method = "long",
    v_ages = v_ages)
  
  expect_equal(summ_prevalence$value, c(rep(0, length(v_ages)-2), prevalence_exp))
  expect_equal(summ_prevalence$person_years_cases, c(rep(0, length(v_ages)-2), n_cohort*prevalence_exp*diff(v_ages[1:2])))
  expect_equal(summ_prevalence$person_years_total, rep(n_cohort*diff(v_ages[1:2]), length(v_ages)-1))
})


###### 4.3 Variation of prevalence
# Set seed for parallelization
set.seed(seed, kind = "L'Ecuyer-CMRG")

# If running locally, use all available cores except for reserved ones
registerDoParallel(cores = detectCores(logical = TRUE) - 2)

stime <- system.time({
  full_summ_prevalence <- foreach(
    i=1:n_sim, 
    .combine=rbind, 
    .inorder=FALSE, 
    .packages=c("data.table","tidyverse")) %dopar% {
      # Simulate cohort
      m_patients <- simulate_data(n_cohort, l_params)
      
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
      
      # Merge longitudinal and cross-sectional
      summ_prevalence <- merge(summ_prevalence,
                               summ_prevalence_long, 
                               by = c("age_start", "age_end"),
                               suffixes = c("_cs", "_long"))
      summ_prevalence
    }
})
print(stime)

# Calculate mean for each age group
mean_prevalence <- full_summ_prevalence[, .(mean_cs = mean(value_cs),
                                            mean_long = mean(value_long)), by = age_start]

# Merge all and mean prevalence
full_summ_prevalence <- merge(full_summ_prevalence,
                              mean_prevalence,
                              by = c("age_start"))

# Check whether CIs contains longitudinal and true prevalence
full_summ_prevalence[, `:=` (contained_true_cs = mean_cs >= ci_lb_cs & mean_cs <= ci_ub_cs,
                             contained_true_long = mean_long >= ci_lb_long & mean_long <= ci_ub_long,
                             consistent_cs_long = value_long >= ci_lb_cs & value_long <= ci_ub_cs)]

# Percentage of values within CIs
pct_contained <- full_summ_prevalence[, .(pct_cs = mean(contained_true_cs),
                                          pct_long = mean(contained_true_long),
                                          pct_consistent = mean(consistent_cs_long)), by = age_start]

# Calculate SD for each age group
sd_prevalence <- full_summ_prevalence[, .(mean_count_cs = mean(n_cases),
                                          mean_total_cs = mean(n_total),
                                          mean_count_long = mean(person_years_cases),
                                          mean_total_long = mean(person_years_total),
                                          mcse_cs = sd(value_cs),
                                          mean_se_cs = mean(se_cs),
                                          mcse_long = sd(value_long),
                                          mean_se_long = mean(se_long)), by = age_start]

# Calculate ratios of SDs
sd_prevalence[, `:=` (ratio_pop = mean_total_long / mean_total_cs,
                      ratio_cs = mcse_cs / mean_se_cs,
                      ratio_long = mcse_long / mean_se_long,
                      ratio_long_cs = mcse_cs / mcse_long)]

# Perform consistency unit tests
test_that("Methods of calculating prevalence produce similar results", {
  expect_equal(abs(pct_contained$pct_cs - conf_level) < 0.03, rep(T, length(v_ages)-1))
  expect_equal(abs(pct_contained$pct_long - conf_level) < 0.03, rep(T, length(v_ages)-1))
})

test_that("Methods of calculating prevalence produce similar results", {
  expect_equal(pct_contained$pct_consistent > conf_level, rep(T, length(v_ages)-1))
})
