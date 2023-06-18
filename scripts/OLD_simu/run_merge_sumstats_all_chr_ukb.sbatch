#!/bin/bash

#SBATCH --job-name=run_merge_sumstats_all_chr_ukb
#SBATCH --cpus-per-task=1
#SBATCH --partition=fm-bigmem-1,fm-bigmem-2,fm-bigmem-3
#SBATCH --time=01:00:00
#SBATCH --mem=15G
#SBATCH --output=run_merge_sumstats_all_chr_ukb.%j.out
#SBATCH --error=run_merge_sumstats_all_chr_ukb.%j.err

###Set some variable
DATA_ID=${1}

###Load module
source /opt/intel/oneapi/mkl/2023.0.0/env/vars.sh intel64
module load R/4.2.3

###Set variables
export MKL_NUM_THREADS=1
export OMP_NUM_THREADS=1

###Run R script
Rscript ../code/merge_sumstats_all_chr_ukb.R --sumstats_prefix ../output/summary_statistics/ukb_caucasian_white_british_unrel_100000_chr \
                                             --sumstats_suffix _sumstats_"${DATA_ID}".rds \
                                             --chr 1:22 \
                                             --output ../output/summary_statistics/ukb_caucasian_white_british_unrel_100000_chrAll_sumstats_"${DATA_ID}".rds \
                                             --seed ${DATA_ID}

###Unload module
module unload R/4.2.3
