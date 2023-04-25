#!/bin/bash

#SBATCH --job-name=compute_LD_sample_filter_ukb_caucasian_white_british_unrel_by_chr
#SBATCH --cpus-per-task=16
#SBATCH --partition=fm-bigmem-1,fm-bigmem-2,fm-bigmem-3
#SBATCH --time=75:00:00
#SBATCH --mem=80G
#SBATCH --output=compute_LD_sample_filter_ukb_caucasian_white_british_unrel_by_chr.%j.out
#SBATCH --error=compute_LD_sample_filter_ukb_caucasian_white_british_unrel_by_chr.%j.err

###Load module
source /opt/intel/oneapi/mkl/2023.0.0/env/vars.sh intel64
module load R/4.2.3

###Set variables
export MKL_NUM_THREADS=1

###Run R script
Rscript ../code/compute_LD_matrix_ukb.R --geno ukb_caucasian_white_british_unrel_100000 \
                                        --test_ids simulated/ukb_caucasian_white_british_unrel_100000_test_ids \
                                        --data_id 1 \
                                        --chr 1:22 \
                                        --impute_missing TRUE \
                                        --ncores 16

###Unload module
module unload R/4.2.3

