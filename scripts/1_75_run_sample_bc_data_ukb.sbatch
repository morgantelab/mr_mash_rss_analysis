#!/bin/bash

#SBATCH --job-name=run_sample_bc_data_ukb
#SBATCH --cpus-per-task=1
#SBATCH --partition=fm-bigmem-1,fm-bigmem-2,fm-bigmem-3
#SBATCH --time=00:30:00
#SBATCH --mem=2G
#SBATCH --output=run_sample_bc_data_ukb.%j.out
#SBATCH --error=run_sample_bc_data_ukb.%j.err

###Load module
source /opt/intel/oneapi/mkl/2023.0.0/env/vars.sh intel64
module load R/4.2.3

###Set variables
export MKL_NUM_THREADS=1
export OMP_NUM_THREADS=1

###Run script
Rscript ../code/sample_bc_data_ukb.R --n_fold 15000

###Unload module
module unload R/4.2.3
