###########################  Vignette for outcome calculations  ###################
#
#  Objective: Vignettes for functions
########################### <<<<<>>>>> ##############################################

rm(list = ls()) # Clean environment
options(scipen = 999) # View data without scientific notation

#### 1.Libraries and functions  ==================================================

###### 1.1 Load packages
library(data.table)
library(dplyr)
library(survival)

###### 1.2 Load functions from R folder
distr.sources <- list.files("R", 
                            pattern="*.R$", full.names=TRUE, 
                            ignore.case=TRUE, recursive = TRUE)
sapply(distr.sources, source, .GlobalEnv)


#### 2. Parameters ========================================================

###### 2.1 General parameters
n_cohort    <- 1000 # Number of individuals
seed <- 123  # Random seed

###### 2.2 Time-to-event parameters
l_params <- list(P  = list(distr = "weibull",
                           params = list(shape = 2,
                                         scale = 300)), # Rate from birth to preclinical cancer onset
                 PC = list(distr = "exp",
                           params = list(rate = 0.2)), # Rate from preclinical to clinical cancer
                 CD = list(distr = "exp",
                           params = list(rate = 0.1)), # Rate from clinical cancer to death
                 Do = list(distr = "weibull",
                           params = list(shape = 5,
                                         scale = 75))) # Rate from birth to death from other causes # Rate from clinical cancer to death

###### 2.3 Epidemiology calculation parameters
method <- "long" # Choose from "cs", "long", or "rcs" for prevalence
v_ages <- seq(40, 80, 10) # Age ranges
rate_unit <- 100000 # Unit for incidence
get_ci <- TRUE # Get confidence intervals


#### 3. Simulate data ========================================================

# Set seed for reproducibility
set.seed(seed)

# Initialize matrix of patient trajectories
m_patients <- cancer_des(n_cohort, l_params)

# Calculate time from cancer to death from any cause
m_patients[, time_CDa := time_D - time_C]


#### 4. Epidemiological calculations ========================================================

# Prevalence
df_prevalence <- calc_prevalence(
  m_patients = m_patients, 
  start_var = "time_P", 
  end_var = "time_C", 
  censor_var = "time_D",
  v_ages = v_ages, 
  method = method,
  output_uncertainty = get_ci)

# Incidence
df_incidence <- calc_incidence(
  m_patients = m_patients, 
  time_var = "time_P", 
  censor_var = "time_D", 
  rate_unit = rate_unit,
  v_ages = v_ages,
  output_uncertainty = get_ci)

# Risk in age range
df_risk <- calc_risk(
  m_patients = m_patients, 
  start_var = "time_P", 
  censor_var = "time_D",
  min_age = min(v_ages), 
  max_age = max(v_ages),
  output_uncertainty = get_ci)

# Mean sojourn time
df_mst <- calc_duration(
  m_patients = m_patients, 
  start_var = "time_P", 
  end_var = "time_C", 
  censor_var = "time_D")

# Survival from diagnosis
df_surv <- calc_surv(
  m_patients = m_patients[time_C < time_D], # Subset of individuals diagnosed before death
  event_time = "time_CDa",
  event_flag = "Dc",
  v_times = seq(0, 10)
)
