#!/bin/bash

#SBATCH --job-name=run_compute_residual_cov_all_chr_ukb_caucasian_white_british_unrel
#SBATCH --cpus-per-task=1
#SBATCH --partition=fm-bigmem-1,fm-bigmem-2,fm-bigmem-3
#SBATCH --time=01:00:00
#SBATCH --mem=4G
#SBATCH --output=run_compute_residual_cov_all_chr_ukb_caucasian_white_british_unrel.%j.out
#SBATCH --error=run_compute_residual_cov_all_chr_ukb_caucasian_white_british_unrel.%j.err

###Set some variable
DATA_ID=${1}

###Load module
source /opt/intel/oneapi/mkl/2023.0.0/env/vars.sh intel64
module load R/4.2.3

###Run R script
Rscript ../code/compute_residual_cov_all_chr_ukb.R --sumstats_prefix ../output/summary_statistics/ukb_caucasian_white_british_unrel_100000_chr \
                                                   --sumstats_suffix _sumstats_"${DATA_ID}".rds \
                                                   --chr 1:22 \
                                                   --output ../output/misc/ukb_caucasian_white_british_unrel_100000_residual_cov_"${DATA_ID}".rds \
                                                   --seed ${DATA_ID}

###Unload module
module unload R/4.2.3

