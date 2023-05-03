#!/bin/bash

#SBATCH --job-name=run_fit_mrmashrss_sample_filter_ukb_caucasian_white_british_unrel_by_chr
#SBATCH --cpus-per-task=4
#SBATCH --partition=fm-bigmem-1,fm-bigmem-2,fm-bigmem-3
#SBATCH --time=12:00:00
#SBATCH --mem=100G
#SBATCH --output=run_fit_mrmashrss_sample_filter_ukb_caucasian_white_british_unrel_by_chr.%j.out
#SBATCH --error=run_fit_mrmashrss_sample_filter_ukb_caucasian_white_british_unrel_by_chr.%j.err

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
DATA_ID=${1}

###Run R script
echo "Start analyzing chromosome ${CHR} ..."

Rscript ../code/fit_mrmashrss_ukb.R --sumstats ukb_caucasian_white_british_unrel_100000_chr"${CHR}"_sumstats \
                                    --LD_matrix ukb_caucasian_white_british_unrel_chr"${CHR}"_LD \
                                    --pheno_cov ukb_caucasian_white_british_unrel_100000_phenotypic_cov \
                                    --residual_cov ukb_caucasian_white_british_unrel_100000_residual_cov \
                                    --data_driven_cov ukb_caucasian_white_british_unrel_100000_bovy_data_driven_cov \
                                    --canonical_cov FALSE \
                                    --n 100000 \
                                    --prop_nonzero 0.01 \
                                    --standardize FALSE \
                                    --update_V_method full \
                                    --w0_threshold 1e-8 \
                                    --X_colmeans ukb_caucasian_white_british_unrel_100000_chr"${CHR}"_geno_means \
                                    --Y_colmeans ukb_caucasian_white_british_unrel_100000_pheno_means \
                                    --ncores 4 \
                                    --data_id ${DATA_ID} \
                                    --outprefix ukb_caucasian_white_british_unrel_100000_chr"${CHR}"

echo "Finished analyzing chromosome ${CHR}."

###Unload module
module unload R/4.2.3
module unload gsl/2.7

