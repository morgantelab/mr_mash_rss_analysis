#!/bin/bash

#SBATCH --job-name=run_compute_prior_all_regions_bc_ukb
#SBATCH --cpus-per-task=10
#SBATCH --partition=fm-bigmem-1,fm-bigmem-2,fm-bigmem-3
#SBATCH --time=24:00:00
#SBATCH --mem=20G
#SBATCH --output=run_compute_prior_all_regions_bc_ukb.%j.out
#SBATCH --error=run_compute_prior_all_regions_bc_ukb.%j.err

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
echo "Start analyzing fold ${FOLD} ..."

Rscript ../code/compute_prior_all_regions_bc_ukb.R --sumstats ../output/summary_statistics/ukb_bc_chrAll_prior_sumstats_"${FOLD}".rds \
                                                   --regions_dir ../data/blood_cell_finemapping_regions \
                                                   --residual_cov full \
                                                   --canonical_cov TRUE \
                                                   --n_strong 1 \
                                                   --n_weak 2 \
                                                   --n_PCs 3 \
                                                   --flash_remove_singleton TRUE \
                                                   --ED_algorithm ${ED} \
                                                   --output ../output/misc/ukb_bc_"${ED}"_prior_"${FOLD}".rds \
                                                   --seed ${FOLD}

echo "Finished analyzing fold ${FOLD}."

###Unload module
module unload R/4.2.3
module unload gsl/2.7
