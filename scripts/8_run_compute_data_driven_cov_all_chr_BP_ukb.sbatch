#!/bin/bash

#SBATCH --job-name=run_compute_data_driven_cov_all_chr_BP_ukb
#SBATCH --cpus-per-task=10
#SBATCH --partition=fm-bigmem-1,fm-bigmem-2,fm-bigmem-3
#SBATCH --time=24:00:00
#SBATCH --mem=20G
#SBATCH --output=run_compute_data_driven_cov_all_chr_BP_ukb.%j.out
#SBATCH --error=run_compute_data_driven_cov_all_chr_BP_ukb.%j.err

###Set some variable
FOLD=${1}
ED=bovy

###Load module
module load gsl/2.7
source /opt/intel/oneapi/mkl/2023.0.0/env/vars.sh intel64
module load R/4.2.3

###Set variables
export MKL_NUM_THREADS=10
export OMP_NUM_THREADS=10

###Run R script
Rscript ../code/compute_data_driven_cov_all_chr_ukb.R --sumstats_prefix ../output/summary_statistics/ukb_tiezzi_BP_chr \
                                                      --sumstats_suffix _sumstats_"${FOLD}".rds \
                                                      --chr 1:22 \
                                                      --residual_cov ../output/misc/ukb_tiezzi_BP_residual_cov_"${FOLD}".rds \
                                                      --strong_Z_thresh 3 \
                                                      --n_PCs 3 \
                                                      --flash_remove_singleton TRUE \
                                                      --ED_algorithm ${ED} \
                                                      --output ../output/misc/ukb_tiezzi_BP_"${ED}"_data_driven_cov_"${FOLD}".rds \
                                                      --seed ${FOLD}

###Unload module
module unload R/4.2.3
module unload gsl/2.7
