#!/bin/bash

#SBATCH --job-name=run_compute_residual_cov_bc_sampled_all_chr_ukb
#SBATCH --cpus-per-task=1
#SBATCH --partition=fm-bigmem-1,fm-bigmem-2,fm-bigmem-3
#SBATCH --time=01:00:00
#SBATCH --mem=4G
#SBATCH --output=run_compute_residual_cov_bc_sampled_all_chr_ukb.%j.out
#SBATCH --error=run_compute_residual_cov_bc_sampled_all_chr_ukb.%j.err

###Set some variable
FOLD=${1}

###Load module
module load gsl/2.7
source /opt/intel/oneapi/mkl/2023.0.0/env/vars.sh intel64
module load R/4.2.3

###Run R script
Rscript ../code/compute_residual_cov_all_chr_ukb.R --sumstats_prefix ../output/summary_statistics/ukb_bc_sampled_chr \
                                                   --sumstats_suffix _sumstats_"${FOLD}".rds \
                                                   --chr 1:22 \
                                                   --output ../output/misc/ukb_bc_sampled_chrAll_residual_cov_"${FOLD}".rds \
                                                   --seed ${FOLD}

###Unload module
module unload R/4.2.3
module load gsl/2.7

