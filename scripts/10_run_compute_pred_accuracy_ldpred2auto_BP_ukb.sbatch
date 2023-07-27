#!/bin/bash

#SBATCH --job-name=run_compute_pred_accuracy_ldpred2auto_BP_ukb
#SBATCH --cpus-per-task=22
#SBATCH --partition=fm-bigmem-1,fm-bigmem-2,fm-bigmem-3
#SBATCH --time=06:00:00
#SBATCH --mem=50G
#SBATCH --output=run_compute_pred_accuracy_ldpred2auto_BP_ukb.%j.out
#SBATCH --error=run_compute_pred_accuracy_ldpred2auto_BP_ukb.%j.err

###Set some variable
FOLD=${1}

###Load module
source /opt/intel/oneapi/mkl/2023.0.0/env/vars.sh intel64
module load R/4.2.3

###Set variables
export MKL_NUM_THREADS=1

###Run R script
Rscript ../code/compute_prediction_accuracy_tiezzi_BP_ukb.R --model ldpred2_auto \
                                                  --model_fit_dir ../output/ldpred2_auto_fit/ \
                                                  --pheno ../data/phenotypes/ukb_tiezzi_cleaned_BP_adjusted_pheno_test_"${FOLD}".rds \
                                                  --geno_dir /data2/morgante_lab/data/ukbiobank/genotypes/imputed/ \
                                                  --sample_file /data2/morgante_lab/data/ukbiobank/genotypes/imputed/ukb22828_c1_b0_v3_s487271.sample \
                                                  --chr 1:22 \
                                                  --trait 1:16 \
                                                  --temp_dir /scratch1/fabiom/tmp \
                                                  --impute_missing FALSE \
                                                  --output_eff ../output/estimated_effects/ukb_tiezzi_BP_ldpred2_auto_effects_"${FOLD}".rds \
                                                  --output_pred_acc ../output/prediction_accuracy/ukb_tiezzi_BP_ldpred2_auto_pred_acc_"${FOLD}".rds \
                                                  --prefix ukb_tiezzi_BP \
                                                  --fold ${FOLD} \
                                                  --ncores 22

###Unload module
module unload R/4.2.3

