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
ED=bovy

###Run R script
echo "Start analyzing chromosome ${CHR} ..."

Rscript ../code/fit_mrmashrss_ukb.R --sumstats ../output/summary_statistics/ukb_caucasian_white_british_unrel_100000_chr"${CHR}"_sumstats_"${DATA_ID}".rds \
                                    --LD_matrix ../data/LD_matrices/ukb_caucasian_white_british_unrel_chr"${CHR}"_LD.ld.bin \
                                    --pheno_cov ../output/misc/ukb_caucasian_white_british_unrel_100000_phenotypic_cov_"${DATA_ID}".rds \
                                    --residual_cov ../output/misc/ukb_caucasian_white_british_unrel_100000_residual_cov_"${DATA_ID}".rds \
                                    --data_driven_cov ../output/misc/ukb_caucasian_white_british_unrel_100000_"${ED}"_data_driven_cov_"${DATA_ID}".rds \
                                    --canonical_cov FALSE \
                                    --n 100000 \
                                    --prop_nonzero 0.01 \
                                    --standardize FALSE \
                                    --update_V_method full \
                                    --w0_threshold 1e-8 \
                                    --X_colmeans ../output/misc/ukb_caucasian_white_british_unrel_100000_chr"${CHR}"_geno_means_"${DATA_ID}".rds \
                                    --Y_colmeans ../output/misc/ukb_caucasian_white_british_unrel_100000_pheno_means_"${DATA_ID}".rds \
                                    --ncores 4 \
                                    --seed ${DATA_ID} \
                                    --outprefix ../output/mr_mash_rss_fit/ukb_caucasian_white_british_unrel_100000_chr"${CHR}"_mr_mash_rss_fit_"${DATA_ID}".rds

echo "Finished analyzing chromosome ${CHR}."

###Unload module
module unload R/4.2.3
module unload gsl/2.7

