#!/bin/bash

#SBATCH --job-name=run_fit_bayesR_by_chr_and_trait_bc_sampled_ukb
#SBATCH --cpus-per-task=1
#SBATCH --partition=fm-bigmem-1,fm-bigmem-2,fm-bigmem-3
#SBATCH --time=24:00:00
#SBATCH --mem=15G
#SBATCH --output=run_fit_bayesR_by_chr_and_trait_bc_sampled_ukb.%j.out
#SBATCH --error=run_fit_bayesR_by_chr_and_trait_bc_sampled_ukb.%j.err

###Load module
source /opt/intel/oneapi/mkl/2023.0.0/env/vars.sh intel64
module load R/4.2.3

###Set variables
export MKL_NUM_THREADS=1
export OMP_NUM_THREADS=1

###Get chromosome number from array task id
CHR="${SLURM_ARRAY_TASK_ID}"

###Set some variables
FOLD=${1}
TRAIT=${2}

###Run R script
echo "Start analyzing trait ${TRAIT}, chromosome ${CHR} and fold ${FOLD} ..."

Rscript ../code/fit_qgg_sbayes_ukb.R --sumstats ../output/summary_statistics/ukb_bc_sampled_chr"${CHR}"_sumstats_"${FOLD}".rds \
                                      --LD_matrix ../data/LD_matrices/ukb_bc_chr"${CHR}"_LD_sparse_"${FOLD}".ld.bin \
                                      --n 60000 \
                                      --method bayesR \
                                      --h2 0.1 \
                                      --pi 0.0001 \
                                      --nburn 1000 \
                                      --nit 5000 \
                                      --nthin 5 \
                                      --trait ${TRAIT} \
                                      --seed ${FOLD} \
                                      --output ../output/bayesR_fit/ukb_bc_sampled_chr"${CHR}"_bayesR_fit_trait"${TRAIT}"_"${FOLD}".rds

echo "Finished analyzing trait ${TRAIT}, chromosome ${CHR} and fold ${FOLD}."

###Unload module
module unload R/4.2.3

