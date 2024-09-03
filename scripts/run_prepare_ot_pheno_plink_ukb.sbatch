#!/bin/bash

#SBATCH --job-name=run_prepare_ot_pheno_plink_ukb
#SBATCH --cpus-per-task=1
#SBATCH --partition=fm-bigmem-1,fm-bigmem-2,fm-bigmem-3,fm-bigmem-4
#SBATCH --time=01:00:00
#SBATCH --mem=4G
#SBATCH --output=run_prepare_ot_pheno_plink_ukb.%j.out
#SBATCH --error=run_prepare_ot_pheno_plink_ukb.%j.err

###Load module
source /opt/intel/oneapi/mkl/2023.0.0/env/vars.sh intel64
module load R/4.2.3

###Set variables
export MKL_NUM_THREADS=1
export OMP_NUM_THREADS=1


Rscript ../code/prepare_ot_pheno_plink_ukb.R --pheno ../data/phenotypes/ukb_cleaned_ot_sampled_covar_pheno_full_data.rds \
                                             --traits "c('waist', 'hip', 'TFM', 'DPa', 'SPa', 'PPa')" \
                                             --output ../data/phenotypes/ukb_cleaned_ot_sampled_adjusted_pheno_full_data_plink.txt

module unload R/4.2.3
