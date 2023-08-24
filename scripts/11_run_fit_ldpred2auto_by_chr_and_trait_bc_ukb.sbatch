#!/bin/bash

#SBATCH --job-name=run_fit_ldpred2auto_by_chr_and_trait_bc_ukb
#SBATCH --cpus-per-task=4
#SBATCH --partition=fm-bigmem-1,fm-bigmem-2,fm-bigmem-3
#SBATCH --time=12:00:00
#SBATCH --mem=350G
#SBATCH --output=run_fit_ldpred2auto_by_chr_and_trait_bc_ukb.%j.out
#SBATCH --error=run_fit_ldpred2auto_by_chr_and_trait_bc_ukb.%j.err

###Load module
source /opt/intel/oneapi/mkl/2023.0.0/env/vars.sh intel64
module load R/4.2.3

###Set variables
export MKL_NUM_THREADS=1

###Get chromosome number from array task id
CHR="${SLURM_ARRAY_TASK_ID}"

###Set some variables
H2_INIT=-1
FOLD=${1}
TRAIT=${2}

###Run R script
echo "Start analyzing trait ${TRAIT}, chromosome ${CHR} and fold ${FOLD} ..."

Rscript ../code/fit_ldpred2auto_ukb.R --sumstats ../output/summary_statistics/ukb_bc_chr"${CHR}"_sumstats_"${FOLD}".rds \
                                      --LD_matrix ../data/LD_matrices/ukb_bc_chr"${CHR}"_LD_sparse.ld.bin \
                                      --n 195239 \
                                      --h2_init "${H2_INIT}" \
                                      --burn_in 500 \
                                      --num_iter 800 \
                                      --trait ${TRAIT} \
                                      --ncores 4 \
                                      --seed ${FOLD} \
                                      --temp_dir /scratch1/fabiom/tmp \
                                      --output ../output/ldpred2_auto_fit/ukb_bc_chr"${CHR}"_ldpred2_auto_fit_trait"${TRAIT}"_"${FOLD}".rds

echo "Finished analyzing trait ${TRAIT}, chromosome ${CHR} and fold ${FOLD}."

###Unload module
module unload R/4.2.3

