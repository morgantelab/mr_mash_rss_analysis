#!/bin/bash

#SBATCH --job-name=run_fit_ldpred2auto_sample_filter_ukb_caucasian_white_british_unrel_gwide_by_trait
#SBATCH --cpus-per-task=4
#SBATCH --partition=fm-bigmem-1,fm-bigmem-2,fm-bigmem-3
#SBATCH --time=24:00:00
#SBATCH --mem=150G
#SBATCH --output=run_fit_ldpred2auto_sample_filter_ukb_caucasian_white_british_unrel_gwide_by_trait.%j.out
#SBATCH --error=run_fit_ldpred2auto_sample_filter_ukb_caucasian_white_british_unrel_gwide_by_trait.%j.err

###Load module
source /opt/intel/oneapi/mkl/2023.0.0/env/vars.sh intel64
module load R/4.2.3

###Set variables
export MKL_NUM_THREADS=1
export OMP_NUM_THREADS=1

###Set some variables
H2_INIT=0.3
DATA_ID=${1}
TRAIT=${2}

###Run R script
echo "Start analyzing trait ${TRAIT} ..."

Rscript ../code/fit_ldpred2auto_gwide_ukb.R --sumstats_prefix ../output/summary_statistics/ukb_caucasian_white_british_unrel_100000_chr \
                                            --sumstats_suffix _sumstats_${DATA_ID}.rds \
                                            --LD_matrix_prefix ../data/LD_matrices/ukb_caucasian_white_british_unrel_chr \
                                            --LD_matrix_suffix _LD_sparse.ld.bin \
                                            --chr 1:22 \
                                            --n 100000 \
                                            --h2_init "${H2_INIT}" \
                                            --burn_in 500 \
                                            --num_iter 500 \
                                            --trait ${TRAIT} \
                                            --ncores 4 \
                                            --seed ${DATA_ID} \
                                            --output ../output/ldpred2_auto_fit/ukb_caucasian_white_british_unrel_100000_chrAll_ldpred2_auto_fit_trait"${TRAIT}"_"${DATA_ID}".rds

echo "Finished analyzing trait ${TRAIT}."

###Unload module
module unload R/4.2.3

