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
n_sim      <- 1000     # Number of simulations

###### 2.2 Time-to-event parameters
l_params_ode <- list(r_P  = 1/200, # Rate from birth to preclinical cancer onset
                     r_PC = 1/10, # Rate from preclinical to clinical cancer
                     r_CD = 1/10, # Rate from clinical cancer to death
                     r_Do = 1/80) # Rate from birth to death from other causes

###### 2.3 Epidemiology calculation parameters
var_onset <- "time_P"
v_ages    <- seq(30, 80, 10) # Age ranges for prevalence


#### 3. Pre-processing  ===========================================

# Set seed for reproducibility
set.seed(seed)


#### 4. Unit tests  ===========================================

###### 4.1 True prevalence - condition ends before death ####
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

# Save column names
v_varnames <- names(m_patients)

# Test accuracy of cross-sectional prevalence
test_that("Test accuracy of cross-sectional prevalence", {
  # Calculate prevalence
  summ_prevalence <- calc_prevalence(
    m_patients, 
    start_var = "time_P", 
    end_var = "time_C", 
    censor_var = "time_D", 
    method = "cs",
    v_ages = v_ages[2:3])
  
  # Check expected value of prevalence
  expect_equal(summ_prevalence$value, prevalence_exp)
  
  # Check expected number of cases
  expect_equal(summ_prevalence$n_cases, n_cohort*prevalence_exp)
  
  # Check expected total population
  expect_equal(summ_prevalence$n_total, n_cohort)
  
  # Check for side effects
  expect_equal(names(m_patients), v_varnames) # Same column names
})

# Test accuracy of longitudinal prevalence
test_that("Test accuracy of longitudinal prevalence", {
  # Calculate prevalence
  summ_prevalence <- calc_prevalence(
    m_patients, 
    start_var = "time_P", 
    end_var = "time_C", 
    censor_var = "time_D", 
    method = "long",
    v_ages = v_ages)
  
  # Check expected value of prevalence
  expect_equal(summ_prevalence$value, rep(prevalence_exp, length(v_ages)-1))
  
  # Check expected person-years with cancer
  expect_equal(summ_prevalence$person_years_cases, rep(n_cohort*prevalence_exp*diff(v_ages[1:2]), length(v_ages)-1))
  
  # Check expected total person-years
  expect_equal(summ_prevalence$person_years_total, rep(n_cohort*diff(v_ages[1:2]), length(v_ages)-1))
})

# Test accuracy of repeated cross-sectional prevalence
test_that("Test accuracy of repeated cross-sectional prevalence", {
  # Calculate prevalence
  summ_prevalence <- calc_prevalence(
    m_patients, 
    start_var = "time_P", 
    end_var = "time_C", 
    censor_var = "time_D", 
    method = "rcs",
    v_ages = v_ages)
  
  # Check expected value of prevalence
  expect_equal(summ_prevalence$value, rep(prevalence_exp, length(v_ages)-1))
  
  # Check expected person-years with cancer
  expect_equal(summ_prevalence$person_years_cases, rep(n_cohort*prevalence_exp*(diff(v_ages[1:2]) + 1), length(v_ages)-1))
  
  # Check expected total person-years
  expect_equal(summ_prevalence$person_years_total, rep(n_cohort*(diff(v_ages[1:2]) + 1), length(v_ages)-1))
})


###### 4.2 True prevalence - condition ends after death ####
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
  
  # Check expected value of prevalence
  expect_equal(summ_prevalence$value, 0)
  
  # Prevalence for last age category should be ~1/4
  summ_prevalence <- calc_prevalence(
    m_patients, 
    start_var = "time_P", 
    end_var = "time_C", 
    censor_var = "time_D", 
    method = "cs",
    v_ages = tail(v_ages, 2))
  
  # Check expected value of prevalence and cases
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
  
  # Check expected value of prevalence, person-years with condition, and person-years total
  expect_equal(summ_prevalence$value, c(rep(0, length(v_ages)-2), prevalence_exp))
  expect_equal(summ_prevalence$person_years_cases, c(rep(0, length(v_ages)-2), n_cohort*prevalence_exp*diff(v_ages[1:2])))
  expect_equal(summ_prevalence$person_years_total, rep(n_cohort*diff(v_ages[1:2]), length(v_ages)-1))
})

# Test accuracy of repeated cross-sectional prevalence
test_that("Test accuracy of repeated cross-sectional prevalence", {
  summ_prevalence <- calc_prevalence(
    m_patients, 
    start_var = "time_P", 
    end_var = "time_C", 
    censor_var = "time_D", 
    method = "rcs",
    v_ages = v_ages)
  
  # Check expected value of prevalence, person-years with condition, and person-years total
  expect_equal(summ_prevalence$value, c(rep(0, length(v_ages)-3), prevalence_exp/(diff(v_ages[1:2])+1), prevalence_exp))
  expect_equal(summ_prevalence$person_years_total, rep(n_cohort*(diff(v_ages[1:2]) + 1), length(v_ages)-1))
})


###### 4.3 Side effects for cross-sectional prevalence ####
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

# Add variables (sample_age and age_start) that should not be used yet in tests
m_patients[, `:=` (sample_age = 0,
                   age_start = 1)]

# Save column names
v_varnames <- names(m_patients)

# Test for cross-sectional prevalence side effects
test_that("Test for cross-sectional prevalence side effects", {
  # Calculate prevalence
  summ_prevalence <- calc_prevalence(
    m_patients, 
    start_var = "time_P", 
    end_var = "time_C", 
    censor_var = "time_D", 
    method = "cs",
    v_ages = v_ages[2:3])
  
  # Check expected value of prevalence
  expect_equal(summ_prevalence$value, prevalence_exp)
  
  # Check expected number of cases
  expect_equal(summ_prevalence$n_cases, n_cohort*prevalence_exp)
  
  # Check expected total population
  expect_equal(summ_prevalence$n_total, n_cohort)
  
  # Check for side effects
  expect_equal(names(m_patients), v_varnames) # Same column names
  expect_equal(m_patients[, sum(sample_age)], 0) # Unchanged sample_age variable
  expect_equal(m_patients[, sum(age_start)], n_cohort) # Unchanged age_start variable
})

# Test that sample age variable works for cross-sectional prevalence
m_patients[, sample_age_var := 65] # Set sample age above eligible age
m_patients[time_P == 0, sample_age_var := mean(v_ages[2:3])] # Set sample age to 45 for people with condition from age 0 to 20
m_patients[m_patients[time_P == v_ages[2], pt_id][1], sample_age_var := mean(v_ages[2:3])] # Set sample age for one individual with condition to 45
var_saved <- m_patients$sample_age_var

# Save column names
v_varnames <- names(m_patients)

test_that("Test sample age for cross-sectional prevalence", {
  # Calculate prevalence
  summ_prevalence <- calc_prevalence(
    m_patients, 
    start_var = "time_P", 
    end_var = "time_C", 
    censor_var = "time_D", 
    method = "cs",
    sample_var = "sample_age_var",
    v_ages = v_ages[2:3])
  
  # Check expected value of prevalence
  expect_equal(summ_prevalence$value, 1 / (n_cohort*prevalence_exp + 1))
  
  # Check expected number of cases
  expect_equal(summ_prevalence$n_cases, 1)
  
  # Check expected total population
  expect_equal(summ_prevalence$n_total, n_cohort*prevalence_exp + 1)
  
  # Check variables (order not necessarily the same)
  expect_true(length(v_varnames[!v_varnames %in% names(m_patients)]) == 0)
  expect_true(length(names(m_patients)[!names(m_patients) %in% v_varnames]) == 0)
  
  # Check that variable values are preserved
  expect_equal(m_patients$sample_age_var, var_saved) # Unchanged sample age variable
  expect_equal(m_patients[, sum(age_start)], n_cohort) # Unchanged age_start variable
})

# Extract sample age as a separate data table, rename sample age variable, and scramble order
dt_sample_ages <- m_patients[, .(pt_id, sample_age = sample_age_var)]
dt_sample_ages <- dt_sample_ages[sample(1:n_cohort, n_cohort, replace = FALSE)]

# Test that (scrambled) sample age data frame works for cross-sectional prevalence
test_that("Test sample age data frame for cross-sectional prevalence", {
  # Calculate prevalence
  summ_prevalence <- calc_prevalence(
    m_patients, 
    start_var = "time_P", 
    end_var = "time_C", 
    censor_var = "time_D", 
    method = "cs",
    dt_sample_ages = dt_sample_ages,
    v_ages = v_ages[2:3])
  
  # Check expected value of prevalence
  expect_equal(summ_prevalence$value, 1 / (n_cohort*prevalence_exp + 1))
  
  # Check expected number of cases
  expect_equal(summ_prevalence$n_cases, 1)
  
  # Check expected total population
  expect_equal(summ_prevalence$n_total, n_cohort*prevalence_exp + 1)
  
  # Check variables (order not necessarily the same)
  expect_true(length(v_varnames[!v_varnames %in% names(m_patients)]) == 0)
  expect_true(length(names(m_patients)[!names(m_patients) %in% v_varnames]) == 0)
  
  # Check that variable values are preserved
  expect_equal(m_patients$sample_age_var, var_saved) # Unchanged sample age variable
  expect_equal(m_patients[, sum(age_start)], n_cohort) # Unchanged age_start variable
})


###### 4.4 Prevalence across all age groups ####
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

# Set sample age variable and save
m_patients[, sample_age := time_P + 1] # Set sample age in range for everyone
m_patients[time_P == v_ages[4], sample_age := max_age + 1] # Reset sample age above max age for fourth group
m_patients[time_P == v_ages[2], sample_age := time_C + 1] # Reset sample age out of disease range for second group
var_saved <- m_patients$sample_age

# Save column names
v_varnames <- names(m_patients)

# Test cross-sectional prevalence across whole population
test_that("Test cross-sectional prevalence across whole population", {
  # Calculate cross-sectional prevalence
  summ_prevalence <- calc_prevalence(
    m_patients, 
    start_var = "time_P", 
    end_var = "time_C", 
    censor_var = "time_D", 
    method = "cs",
    sample_var = "sample_age") # Note: v_ages is NULL
  
  # Check expected value of prevalence
  expect_equal(summ_prevalence$value, (1/prevalence_exp-2)/(1/prevalence_exp-1))
  
  # Check expected number of cases
  expect_equal(summ_prevalence$n_cases, n_cohort*(1/prevalence_exp-2)*prevalence_exp)
  
  # Check expected total population
  expect_equal(summ_prevalence$n_total, n_cohort*(1-prevalence_exp))
  
  # Check for side effects
  expect_equal(names(m_patients), v_varnames) # Same column names
  expect_equal(m_patients$sample_age, var_saved) # Unchanged sample age variable
})

# Test longitudinal prevalence across whole population
test_that("Test longitudinal prevalence across whole population", {
  # Calculate longitudinal prevalence
  summ_prevalence <- calc_prevalence(
    m_patients, 
    start_var = "time_P", 
    end_var = "time_C", 
    censor_var = "time_D", 
    method = "long") # Note: v_ages is NULL
  
  # Check expected value of prevalence
  expect_equal(summ_prevalence$value, prevalence_exp)
  
  # Check expected person-years with cancer
  expect_equal(summ_prevalence$person_years_cases, n_cohort*prevalence_exp*max_age)
  
  # Check expected total person-years
  expect_equal(summ_prevalence$person_years_total, n_cohort*max_age)
})


###### 4.5 Variation of prevalence
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

# Map age range groups
df_cancer_cohort_wide[, `:=` (age_idx = findInterval(time, v_ages),
                              age_idx_ub = findInterval(time, v_ages, left.open = T))]
df_cancer_cohort_wide[, `:=` (age_start = v_ages[age_idx]), by = age_idx]
df_cancer_cohort_wide[age_idx != age_idx_ub, `:=` (age_grp_ub = v_ages[age_idx_ub]), by = age_idx_ub]

# Solve for true prevalence in each age range
true_prevalence <- df_cancer_cohort_wide[age_start < max(v_ages), 
                                         .(p_cases = sum(P),
                                           p_total = sum(H, P, C)),
                                         by = age_start]

# Extract boundary counts
true_prevalence_ub <- df_cancer_cohort_wide[!is.na(age_grp_ub)][, p_total := H + P + C]

# Add boundary counts to each age group
true_prevalence[, `:=` (p_cases = p_cases + true_prevalence_ub$P,
                        p_total = p_total + true_prevalence_ub$p_total)]

# Calculate true prevalence
true_prevalence[, `:=` (true_value = p_cases / p_total)]
v_true <- true_prevalence$true_value

# Run simulations
stime <- system.time({
  full_summ_prevalence <- foreach(
    i=1:n_sim, 
    .combine=rbind, 
    .inorder=FALSE, 
    .packages=c("data.table","tidyverse")) %dopar% {
      # Simulate cohort
      m_patients <- cancer_des(n_cohort, l_params_des)
      
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
full_summ_prevalence[, `:=` (
  contained_true_cs = v_true >= ci_lb & v_true <= ci_ub,
  # contained_true_long = v_true >= ci_lb_long & v_true <= ci_ub_long,
  # contained_true_rcs = v_true >= ci_lb & v_true <= ci_ub,
  consistent_cs_long = value_long >= ci_lb & value_long <= ci_ub,
  consistent_cs_rcs = value >= ci_lb & value <= ci_ub
)]

# Percentage of values within CIs
pct_contained <- full_summ_prevalence[, .(pct_cs = mean(contained_true_cs),
                                          # pct_long = mean(contained_true_long),
                                          # pct_rcs = mean(contained_true_rcs),
                                          pct_consistent_long = mean(consistent_cs_long),
                                          pct_consistent_rcs = mean(consistent_cs_rcs)), 
                                      by = age_start]

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
