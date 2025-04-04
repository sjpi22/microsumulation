###########################  Vignette for outcome calculations  ###################
#
#  Objective: Vignettes for functions
########################### <<<<<>>>>> ##############################################

rm(list = ls()) # Clean environment
options(scipen = 999) # View data without scientific notation

#### 1.Libraries and functions  ==================================================

###### 1.1 Load packages
library(data.table)

###### 1.2 Load functions from R folder
distr.sources <- list.files("R", 
                            pattern="*.R$", full.names=TRUE, 
                            ignore.case=TRUE, recursive = TRUE)
sapply(distr.sources, source, .GlobalEnv)


#### 2. Parameters ========================================================

###### 2.1 General parameters
n    <- 1000 # Number of individuals
seed <- 123  # Random seed

###### 2.2 Time-to-event parameters
params_H_Do <- c(min = 10, max = 50)  # Time from birth to death from other causes
params_H_S  <- c(rate = 0.1)          # Time from birth to disease onset
params_S_R  <- c(rate = 0.3)            # Time from disease onset to recovery
params_S_Dx <- c(rate = 0.1)            # Time from disease onset to death from disease

###### 2.3 Epidemiology calculation parameters
v_ages <- seq(0, 50, 10) # Age ranges for prevalence
  
  
#### 3. Simulate data ========================================================

# Set seed for reproducibility
set.seed(seed)

# Initialize matrix of patient trajectories
m_patients <- data.table(pt_id = 1:n)
setkey(m_patients, pt_id)

# Simulate time to disease onset
m_patients[, `:=` (time_H_Do = runif(.N, min = params_H_Do["min"], max = params_H_Do["max"]),
                   time_H_S  = rexp(.N, rate = params_H_S["rate"]),
                   time_S_R  = rexp(.N, rate = params_S_R["rate"]),
                   time_S_Dx = rexp(.N, rate = params_S_Dx["rate"]))]

# Calculate time from birth to recovery and death from disease
m_patients[time_S_R < time_S_Dx, `:=` (time_H_R = time_H_S + time_S_R)]
m_patients[time_S_R >= time_S_Dx, `:=` (time_H_Dx = time_H_S + time_S_Dx)]

# Calculate time to death from any cause
m_patients[, time_H_D := pmin(time_H_Do, time_H_Dx, na.rm = T)]

# Calculate disease start and end time if disease occurred during lifetime
m_patients[time_H_S < time_H_D, `:=` (time_start = time_H_S,
                                      time_end   = pmin(time_H_R, time_H_D, na.rm = T))]


#### 4. Epidemiological calculations ========================================================

df_prev <- calc_prevalence(m_patients = m_patients, 
                           start_var = "time_start", 
                           end_var = "time_end", 
                           censor_var = "time_H_D", 
                           v_ages = v_ages)

