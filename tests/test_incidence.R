###########################  Unit test: Incidence   ################
#
#  Objective: Run unit tests for calculating incidence
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
    method = "long",
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
    method = "long",
    rate_unit = rate_unit)
  
  # Check expected value of incidence
  expect_equal(summ_incidence$value, c(0, 1/22)*rate_unit)
  
  # Check expected person-years with cancer
  expect_equal(summ_incidence$n_events, c(0, 1))
  
  # Check expected total person-years
  expect_equal(summ_incidence$person_years_total, c(40, 22))
})


###### 4.2 Variation of incidence
# Set seed for parallelization
set.seed(seed, kind = "L'Ecuyer-CMRG")

# If running locally, use all available cores except for reserved ones
registerDoParallel(cores = detectCores(logical = TRUE) - 2)

stime <- system.time({
  full_summ_incidence <- foreach(
    i=1:n_sim, 
    .combine=rbind, 
    .inorder=FALSE, 
    .packages=c("data.table","tidyverse")) %dopar% {
      # Simulate cohort
      m_patients <- cancer_des(n_cohort, l_params)
      
      # Calculate incidence (longitudinal)
      summ_incidence_long <- calc_incidence(
        m_patients, 
        time_var = "time_C", 
        censor_var = "time_D", 
        method = "long",
        v_ages = v_ages,
        output_uncertainty = T,
        rate_unit = rate_unit)
      
      summ_incidence_long
    }
})
print(stime)

# Calculate mean for each age group
mean_incidence <- full_summ_incidence[, .(mean_long = mean(value)), by = age_start]

# Merge all and mean incidence
full_summ_incidence <- merge(full_summ_incidence,
                              mean_incidence,
                              by = c("age_start"))

# Check whether CIs contains longitudinal and true incidence
full_summ_incidence[, `:=` (contained_true_long = mean_long >= ci_lb & mean_long <= ci_ub)]

# Percentage of values within CIs
pct_contained <- full_summ_incidence[, .(pct_long = mean(contained_true_long)), by = age_start]

# Calculate SD for each age group
sd_incidence <- full_summ_incidence[, .(mcse_long = sd(value),
                                        mean_se_long = mean(se)), by = age_start]

# Calculate ratios of SDs
sd_incidence[, `:=` (ratio_long = mcse_long / mean_se_long)]

# Perform consistency unit tests
test_that("Methods of calculating incidence produce similar results", {
  expect_equal(abs(pct_contained$pct_long - conf_level) < 0.03, rep(T, length(v_ages)-1))
})

