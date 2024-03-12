#!/bin/bash

#SBATCH --job-name=run_ukb_bc_all_data_plot
#SBATCH --cpus-per-task=1
#SBATCH --partition=fm-bigmem-1,fm-bigmem-2,fm-bigmem-3
#SBATCH --time=00:30:00
#SBATCH --mem=2G
#SBATCH --output=run_ukb_bc_all_data_plot.%j.out
#SBATCH --error=run_ukb_bc_all_data_plot.%j.err

###Load module
source /opt/intel/oneapi/mkl/2023.0.0/env/vars.sh intel64
module load R/4.2.3

###Set variables
export MKL_NUM_THREADS=1
export OMP_NUM_THREADS=1

###Run script
Rscript ../code/paper_figures/ukb_bc_all_data_plot.R

###Unload module
module unload R/4.2.3
