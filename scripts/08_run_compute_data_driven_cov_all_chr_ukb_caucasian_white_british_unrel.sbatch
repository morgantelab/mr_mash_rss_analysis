#!/bin/bash

#SBATCH --job-name=run_compute_data_driven_cov_all_chr_ukb_caucasian_white_british_unrel
#SBATCH --cpus-per-task=8
#SBATCH --partition=fm-bigmem-1,fm-bigmem-2,fm-bigmem-3
#SBATCH --time=24:00:00
#SBATCH --mem=7G
#SBATCH --output=run_compute_data_driven_cov_all_chr_ukb_caucasian_white_british_unrel.%j.out
#SBATCH --error=run_compute_data_driven_cov_all_chr_ukb_caucasian_white_british_unrel.%j.err

###Set some variable
DATA_ID=${1}

###Load module
module load gsl/2.7
source /opt/intel/oneapi/mkl/2023.0.0/env/vars.sh intel64
module load R/4.2.3

###Set variables
export MKL_NUM_THREADS=8
export OMP_NUM_THREADS=8

###Run R script
Rscript ../code/compute_data_driven_cov_all_chr_ukb.R --sumstats ukb_caucasian_white_british_unrel_100000 \
                                                      --residual_cov ukb_caucasian_white_british_unrel_100000_residual_cov \
                                                      --strong_Z_thresh 3 \
                                                      --n_PCs 3 \
                                                      --flash_remove_singleton TRUE \
                                                      --ED_algorithm bovy \
                                                      --data_id ${DATA_ID}

###Unload module
module unload R/4.2.3
module unload gsl/2.7
