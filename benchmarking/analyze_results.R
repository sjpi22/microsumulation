###########################  Benchmarking analysis  ################
#
#  Objective: Run benchmarking simulations
########################### <<<<<>>>>> ##############################################

rm(list = ls()) # Clean environment
options(scipen = 999) # View data without scientific notation

#### 1.Libraries and functions  ==================================================

###### 1.1 Load packages
library(tidyverse)
library(data.table)

###### 1.2 Load functions

# Load functions
distr.sources <- list.files("R", 
                            pattern="*.R$", full.names=TRUE, 
                            ignore.case=TRUE, recursive = TRUE)
sapply(distr.sources, source, .GlobalEnv)


#### 2. General parameters ========================================================

###### 2.1 General parameters
path_benchmark <- "benchmarking/result.rds" # Path to save results

###### 2.3 Epidemiology calculation parameters
l_age_exp <- list(
  c(30, 80),
  seq(30, 80, 10) # Age ranges for prevalence
)


#### 3. Pre-processing  ===========================================

# Load results
res <- readRDS(path_benchmark)


#### 4. Analysis  ===========================================


res$time[[2]]
test = res$outcomes[[2]]
test[, sd_ratio := sd / mean]
