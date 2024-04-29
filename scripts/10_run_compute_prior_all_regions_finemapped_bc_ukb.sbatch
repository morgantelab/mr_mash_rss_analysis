#!/bin/bash

#SBATCH --job-name=run_compute_prior_all_regions_finemapped_bc_ukb
#SBATCH --cpus-per-task=16
#SBATCH --partition=fm-bigmem-4,fm-bigmem-1,fm-bigmem-2,fm-bigmem-3
#SBATCH --time=96:00:00
#SBATCH --mem=200G
#SBATCH --output=run_compute_prior_all_regions_finemapped_bc_ukb.%j.out
#SBATCH --error=run_compute_prior_all_regions_finemapped_bc_ukb.%j.err

###Set some variable
FOLD=${1}
ED=bovy

###Load module
module load gsl/2.7
source /opt/intel/oneapi/mkl/2023.0.0/env/vars.sh intel64
module load R/4.2.3

###Set variables
export MKL_NUM_THREADS=16
export OMP_NUM_THREADS=16

###Run R script
echo "Start analyzing fold ${FOLD} ..."

Rscript ../code/compute_prior_all_regions_finemapped_bc_ukb.R --sumstats ../output/summary_statistics/ukb_bc_chrAll_prior_sumstats_"${FOLD}".rds \
                                                              --regions_dir ../data/blood_cell_finemapping_regions \
                                                              --residual_cov full \
                                                              --canonical_cov TRUE \
                                                              --geno /scratch3/fabiom/ukbiobank_files/ukb_bc_geno_imp_prior.rds \
                                                              --fold_file ../data/phenotypes/ukb_cleaned_bc_covar_pheno.rds \
                                                              --sample_file /data2/morgante_lab/data/ukbiobank/genotypes/imputed/ukb22828_c1_b0_v3_s487271.sample \
                                                              --n_weak 2 \
                                                              --n_PCs 3 \
                                                              --flash_remove_singleton TRUE \
                                                              --ED_algorithm ${ED} \
                                                              --fold ${FOLD} \
                                                              --output ../output/misc/ukb_bc_"${ED}"_prior_finemapped_"${FOLD}".rds \
                                                              --seed ${FOLD} \
                                                             --temp_dir /scratch3/fabiom/tmp

echo "Finished analyzing fold ${FOLD}."

###Unload module
module unload R/4.2.3
module unload gsl/2.7
