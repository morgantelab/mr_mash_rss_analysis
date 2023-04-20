#!/bin/bash

#SBATCH --job-name=run_compute_compute_residual_cov_all_chr_ukb_caucasian_white_british_unrel
#SBATCH --cpus-per-task=1
#SBATCH --partition=fm-bigmem-1,fm-bigmem-2,fm-bigmem-3
#SBATCH --time=01:00:00
#SBATCH --mem=4G
#SBATCH --output=run_compute_compute_residual_cov_all_chr_ukb_caucasian_white_british_unrel.%j.out
#SBATCH --error=run_compute_compute_residual_cov_all_chr_ukb_caucasian_white_british_unrel.%j.err

###Load module
source /opt/intel/oneapi/mkl/2023.0.0/env/vars.sh intel64
module load R/4.2.3

###Run R script
Rscript ../code/compute_residual_cov_all_chr_ukb.R --input ukb_caucasian_white_british_unrel_100000 \
                                       --data_id 1

###Unload module
module unload R/4.2.3

