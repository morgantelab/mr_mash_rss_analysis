#!/bin/bash

#SBATCH --job-name=run_fit_bayesR_sample_filter_ukb_caucasian_white_british_unrel_by_chr_and_trait
#SBATCH --cpus-per-task=4
#SBATCH --partition=fm-bigmem-1,fm-bigmem-2,fm-bigmem-3
#SBATCH --time=24:00:00
#SBATCH --mem=100G
#SBATCH --output=run_fit_bayesR_sample_filter_ukb_caucasian_white_british_unrel_by_chr_and_trait.%j.out
#SBATCH --error=run_fit_bayesR_sample_filter_ukb_caucasian_white_british_unrel_by_chr_and_trait.%j.err

###Load module
source /opt/intel/oneapi/mkl/2023.0.0/env/vars.sh intel64
module load R/4.2.3

###Set variables
export MKL_NUM_THREADS=4
export OMP_NUM_THREADS=4

###Get chromosome number from array task id
CHR="${SLURM_ARRAY_TASK_ID}"

###Set some variables
H2_INIT=0.013
DATA_ID=${1}
TRAIT=${2}

###Run R script
echo "Start analyzing trait ${TRAIT} and chromosome ${CHR} ..."

Rscript ../code/fit_qgg_sbayes_ukb.R --sumstats ../output/summary_statistics/ukb_caucasian_white_british_unrel_100000_chr"${CHR}"_sumstats_"${DATA_ID}".rds \
                                      --LD_matrix ../data/LD_matrices/ukb_caucasian_white_british_unrel_chr"${CHR}"_LD.ld.bin \
                                      --n 100000 \
                                      --method bayesR \
                                      --h2 "${H2_INIT}" \
                                      --nburn 200 \
                                      --nit 1000 \
                                      --updateB TRUE \
                                      --updateE TRUE \
                                      --updatePi TRUE \
                                      --updateG TRUE \
                                      --adjustE FALSE \
                                      --trait ${TRAIT} \
                                      --seed ${DATA_ID} \
                                      --output ../output/bayesR_fit/ukb_caucasian_white_british_unrel_100000_chr"${CHR}"_bayesR_fit_trait"${TRAIT}"_"${DATA_ID}".rds

echo "Finished analyzing trait ${TRAIT} and chromosome ${CHR}."

###Unload module
module unload R/4.2.3

