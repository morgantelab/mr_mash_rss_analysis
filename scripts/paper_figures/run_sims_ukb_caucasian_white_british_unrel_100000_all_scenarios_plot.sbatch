#!/bin/bash

#SBATCH --job-name=run_sims_ukb_caucasian_white_british_unrel_100000_all_scenarios_plot
#SBATCH --cpus-per-task=1
#SBATCH --partition=fm-bigmem-1,fm-bigmem-2,fm-bigmem-3
#SBATCH --time=00:30:00
#SBATCH --mem=2G
#SBATCH --output=run_sims_ukb_caucasian_white_british_unrel_100000_all_scenarios_plot.%j.out
#SBATCH --error=run_sims_ukb_caucasian_white_british_unrel_100000_all_scenarios_plot.%j.err

###Load module
source /opt/intel/oneapi/mkl/2023.0.0/env/vars.sh intel64
module load R/4.2.3

###Set variables
export MKL_NUM_THREADS=1
export OMP_NUM_THREADS=1

###Run script
Rscript ../code/paper_figures/sims_ukb_caucasian_white_british_unrel_100000_all_scenarios_plot.R

###Unload module
module unload R/4.2.3
