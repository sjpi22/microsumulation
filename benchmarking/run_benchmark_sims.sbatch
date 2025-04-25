#!/usr/bin/bash
#SBATCH --job-name=run_benchmark_sims
#SBATCH --nodes=1
#SBATCH --time=2:00:00
#SBATCH --mem=3G
#SBATCH --output=benchmarking/run_benchmark_sims.log
#SBATCH --mail-user=sjpi@stanford.edu
#SBATCH --mail-type=ALL

#* To send the job to Sherlock use:
# - sbatch benchmarking/run_benchmark_sims.sbatch

# load R module
ml R/4.3.2

# run R code
Rscript benchmarking/run_benchmark_sims.R
