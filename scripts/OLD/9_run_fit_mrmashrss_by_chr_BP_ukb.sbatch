#!/bin/bash

#SBATCH --job-name=run_fit_mrmashrss_by_chr_BP_ukb
#SBATCH --cpus-per-task=4
#SBATCH --partition=fm-bigmem-1,fm-bigmem-2,fm-bigmem-3
#SBATCH --time=72:00:00
#SBATCH --mem=200G
#SBATCH --output=run_fit_mrmashrss_by_chr_BP_ukb.%j.out
#SBATCH --error=run_fit_mrmashrss_by_chr_BP_ukb.%j.err

###Load module
module load gsl/2.7
source /opt/intel/oneapi/mkl/2023.0.0/env/vars.sh intel64
module load R/4.2.3

###Set variables
export MKL_NUM_THREADS=1
export OMP_NUM_THREADS=1

###Get chromosome number from array task id
CHR="${SLURM_ARRAY_TASK_ID}"

###Set some variable
FOLD=${1}
ED=bovy

###Run R script
echo "Start analyzing chromosome ${CHR} and fold ${FOLD} ..."

Rscript ../code/fit_mrmashrss_ukb.R --sumstats ../output/summary_statistics/ukb_tiezzi_BP_chr"${CHR}"_sumstats_"${FOLD}".rds \
                                    --LD_matrix ../data/LD_matrices/ukb_caucasian_white_british_unrel_for_tiezzi_BP_chr"${CHR}"_LD.ld.bin \
                                    --pheno_cov ../output/misc/ukb_tiezzi_BP_phenotypic_cov_"${FOLD}".rds \
                                    --residual_cov ../output/misc/ukb_tiezzi_BP_residual_cov_"${FOLD}".rds \
                                    --data_driven_cov ../output/misc/ukb_tiezzi_BP_"${ED}"_data_driven_cov_"${FOLD}".rds \
                                    --canonical_cov TRUE \
                                    --n 75311 \
                                    --prop_nonzero 0.01 \
                                    --standardize FALSE \
                                    --update_V_method full \
                                    --w0_threshold 1e-8 \
                                    --ncores 4 \
                                    --seed ${FOLD} \
                                    --output ../output/mr_mash_rss_fit/ukb_tiezzi_BP_chr"${CHR}"_mr_mash_rss_fit_"${FOLD}".rds

echo "Finished analyzing chromosome ${CHR} and fold ${FOLD}."

###Unload module
module unload R/4.2.3
module unload gsl/2.7

