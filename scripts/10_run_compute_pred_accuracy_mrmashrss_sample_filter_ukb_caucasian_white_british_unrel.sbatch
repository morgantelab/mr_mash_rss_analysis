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
                                                  --model_fit_dir ../output/mr_mash_rss_fit/ \
                                                  --pheno ../data/phenotypes/simulated/ukb_caucasian_white_british_unrel_100000_pheno_"${DATA_ID}".rds \
                                                  --geno ../data/genotypes/ukb_caucasian_white_british_unrel_100000.bed \
                                                  --pheno_means ../output/misc/ukb_caucasian_white_british_unrel_100000_pheno_means_"${DATA_ID}".rds \
                                                  --test_ids ../data/phenotypes/simulated/ukb_caucasian_white_british_unrel_100000_test_ids_"${DATA_ID}".txt \
                                                  --data_id ${DATA_ID} \
                                                  --chr 1:22 \
                                                  --trait 1:5 \
                                                  --temp_dir /data2/morgante_lab/fabiom/tmp \
                                                  --impute_missing TRUE \
                                                  --output_eff ../output/estimated_effects/ukb_caucasian_white_british_unrel_100000_mr_mash_rss_effects_"${DATA_ID}".rds \
                                                  --output_pred_acc ../output/prediction_accuracy/ukb_caucasian_white_british_unrel_100000_mr_mash_rss_pred_acc_"${DATA_ID}".rds \
                                                  --ncores 4

###Unload module
module unload R/4.2.3

