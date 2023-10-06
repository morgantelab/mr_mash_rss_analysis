#!/bin/bash

#SBATCH --job-name=run_fit_ldsc_gwide_by_trait_bc_uk
#SBATCH --cpus-per-task=4
#SBATCH --partition=fm-bigmem-1,fm-bigmem-2,fm-bigmem-3
#SBATCH --time=24:00:00
#SBATCH --mem=400G
#SBATCH --output=run_fit_ldsc_gwide_by_trait_bc_uk.%j.out
#SBATCH --error=run_fit_ldsc_gwide_by_trait_bc_uk.%j.err
#SBATCH --array=1-16

###Load module
source /opt/intel/oneapi/mkl/2023.0.0/env/vars.sh intel64
module load R/4.2.3

###Set variables
export MKL_NUM_THREADS=1

###Set some variables
FOLD=${1}
TRAIT="${SLURM_ARRAY_TASK_ID}"

###Run R script
echo "Start analyzing trait ${TRAIT} and fold ${FOLD} ..."

Rscript ../code/fit_ldsc_ukb.R --sumstats_prefix ../output/summary_statistics/ukb_bc_chr \
                               --sumstats_suffix _sumstats_"${FOLD}".rds \
                               --LD_matrix_prefix ../data/LD_matrices/ukb_bc_chr \
                               --LD_matrix_suffix _LD_sparse_"${FOLD}".ld.bin \
                               --chr 1:22 \
                               --n 195239 \
                               --trait ${TRAIT} \
                               --ncores 4 \
                               --seed ${FOLD} \
                               --output ../output/ldsc_fit/ukb_bc_chrAll_ldsc_fit_trait"${TRAIT}"_"${FOLD}".rds

echo "Finished analyzing trait ${TRAIT} and fold ${FOLD}."

###Unload module
module unload R/4.2.3

