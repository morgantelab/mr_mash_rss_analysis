#!/bin/bash

#SBATCH --job-name=run_compute_pred_accuracy_mrmashrss_sample_filter_ukb_caucasian_white_british_unrel
#SBATCH --cpus-per-task=4
#SBATCH --partition=fm-bigmem-1,fm-bigmem-2,fm-bigmem-3
#SBATCH --time=02:00:00
#SBATCH --mem=100G
#SBATCH --output=run_compute_pred_accuracy_mrmashrss_sample_filter_ukb_caucasian_white_british_unrel.%j.out
#SBATCH --error=run_compute_pred_accuracy_mrmashrss_sample_filter_ukb_caucasian_white_british_unrel.%j.err

###Set some variable
DATA_ID=${1}

###Load module
source /opt/intel/oneapi/mkl/2023.0.0/env/vars.sh intel64
module load R/4.2.3

###Set variables
export MKL_NUM_THREADS=1

###Run R script
Rscript ../code/compute_prediction_accuracy_ukb.R --model mr_mash_rss \
                                                  --pheno simulated/ukb_caucasian_white_british_unrel_100000_pheno \
                                                  --geno ukb_caucasian_white_british_unrel_100000 \
                                                  --pheno_means ukb_caucasian_white_british_unrel_100000_pheno_means \
                                                  --test_ids simulated/ukb_caucasian_white_british_unrel_100000_test_ids \
                                                  --prefix ukb_caucasian_white_british_unrel_100000 \
                                                  --data_id ${DATA_ID} \
                                                  --chr 1:22 \
                                                  --trait 1:5 \
                                                  --impute_missing TRUE \
                                                  --ncores 4

###Unload module
module unload R/4.2.3

