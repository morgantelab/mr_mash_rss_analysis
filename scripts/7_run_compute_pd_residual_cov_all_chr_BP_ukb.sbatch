#!/bin/bash

#SBATCH --job-name=run_compute_pd_residual_cov_all_chr_BP_ukb
#SBATCH --cpus-per-task=1
#SBATCH --partition=fm-bigmem-1,fm-bigmem-2,fm-bigmem-3
#SBATCH --time=01:00:00
#SBATCH --mem=4G
#SBATCH --output=run_compute_pd_residual_cov_all_chr_BP_ukb.%j.out
#SBATCH --error=run_compute_pd_residual_cov_all_chr_BP_ukb.%j.err

###Set some variable
FOLD=${1}

###Load module
source /opt/intel/oneapi/mkl/2023.0.0/env/vars.sh intel64
module load R/4.2.3

###Run R script
Rscript ../code/compute_pd_residual_cov_all_chr_ukb.R --sumstats_prefix ../output/summary_statistics/ukb_tiezzi_BP_chr \
                                                      --sumstats_suffix _sumstats_"${FOLD}".rds \
                                                      --chr 1:22 \
                                                      --output ../output/misc/ukb_tiezzi_BP_pd_residual_cov_"${FOLD}".rds \
                                                      --seed ${FOLD}

###Unload module
module unload R/4.2.3

