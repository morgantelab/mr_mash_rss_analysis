#!/bin/bash

#SBATCH --job-name=run_compute_pred_accuracy_mrmashrss_sparse_LD_V_all_chr_init_prior_finemapped_bc_sampled_ukb
#SBATCH --cpus-per-task=8
#SBATCH --partition=fm-bigmem-1,fm-bigmem-2,fm-bigmem-3
#SBATCH --time=06:00:00
#SBATCH --mem=25G
#SBATCH --output=run_compute_pred_accuracy_mrmashrss_sparse_LD_V_all_chr_init_prior_finemapped_bc_sampled_ukb.%j.out
#SBATCH --error=run_compute_pred_accuracy_mrmashrss_sparse_LD_V_all_chr_init_prior_finemapped_bc_sampled_ukb.%j.err

###Set some variable
FOLD=${1}

###Load module
source /opt/intel/oneapi/mkl/2023.0.0/env/vars.sh intel64
module load R/4.2.3

###Set variables
export MKL_NUM_THREADS=1

###Run R script
Rscript ../code/compute_prediction_accuracy_bc_ukb.R --model mr_mash_rss_sparse_LD_V_all_chr_init_prior_finemapped \
                                                     --model_fit_dir ../output/mr_mash_rss_fit/ \
                                                     --pheno ../data/phenotypes/ukb_cleaned_bc_sampled_adjusted_pheno_test_"${FOLD}".rds \
                                                     --geno /scratch1/fabiom/ukb_bc_sampled_geno_imp_HM3.rds \
                                                     --sample_file /data2/morgante_lab/data/ukbiobank/genotypes/imputed/ukb22828_c1_b0_v3_s487271.sample \
                                                     --samples_to_keep ../data/misc/ukb_cleaned_bc_sampled_ind_ids.txt \
                                                     --chr 1:22 \
                                                     --trait 1:16 \
                                                     --impute_missing FALSE \
                                                     --output_eff ../output/estimated_effects/ukb_bc_sampled_mr_mash_rss_sparse_LD_V_all_chr_init_prior_finemapped_effects_"${FOLD}".rds \
                                                     --output_pred_acc ../output/prediction_accuracy/ukb_bc_sampled_mr_mash_rss_sparse_LD_V_all_chr_init_prior_finemapped_pred_acc_"${FOLD}".rds \
                                                     --prefix ukb_bc_sampled \
                                                     --fold ${FOLD} \
                                                     --ncores 8

###Unload module
module unload R/4.2.3

