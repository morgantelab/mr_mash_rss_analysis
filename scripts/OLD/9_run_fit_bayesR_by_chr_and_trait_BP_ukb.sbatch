#!/bin/bash

#SBATCH --job-name=run_fit_bayesR_by_chr_and_trait_BP_ukb
#SBATCH --cpus-per-task=1
#SBATCH --partition=bigmem
#SBATCH --time=24:00:00
#SBATCH --mem=150G
#SBATCH --output=run_fit_bayesR_by_chr_and_trait_BP_ukb.%j.out
#SBATCH --error=run_fit_bayesR_by_chr_and_trait_BP_ukb.%j.err

###Load module
source /opt/intel/oneapi/mkl/2023.0.0/env/vars.sh intel64
module load R/4.2.3

###Set variables
export MKL_NUM_THREADS=1
export OMP_NUM_THREADS=1

###Get chromosome number from array task id
CHR="${SLURM_ARRAY_TASK_ID}"

###Set some variables
H2_INIT=0.009
FOLD=${1}
TRAIT=${2}

###Run R script
echo "Start analyzing trait ${TRAIT}, chromosome ${CHR} and fold ${FOLD} ..."

Rscript ../code/fit_qgg_sbayes_ukb.R --sumstats ../output/summary_statistics/ukb_tiezzi_BP_chr"${CHR}"_sumstats_"${FOLD}".rds \
                                      --LD_matrix ../data/LD_matrices/ukb_caucasian_white_british_unrel_for_tiezzi_BP_chr"${CHR}"_LD.ld.bin \
                                      --n 75311 \
                                      --method bayesR \
                                      --h2 "${H2_INIT}" \
                                      --nburn 200 \
                                      --nit 800 \
                                      --nub 4 \
                                      --nue 4 \
                                      --updateB TRUE \
                                      --updateE TRUE \
                                      --updatePi TRUE \
                                      --updateG TRUE \
                                      --adjustE TRUE \
                                      --trait ${TRAIT} \
                                      --seed ${FOLD} \
                                      --output ../output/bayesR_fit/ukb_tiezzi_BP_chr"${CHR}"_bayesR_fit_trait"${TRAIT}"_"${FOLD}".rds

echo "Finished analyzing trait ${TRAIT}, chromosome ${CHR} and fold ${FOLD}."

###Unload module
module unload R/4.2.3

