#!/bin/bash

#SBATCH --job-name=run_fit_ldpred2auto_sample_filter_ukb_caucasian_white_british_unrel_by_chr_and_trait
#SBATCH --cpus-per-task=4
#SBATCH --partition=fm-bigmem-1,fm-bigmem-2,fm-bigmem-3
#SBATCH --time=12:00:00
#SBATCH --mem=100G
#SBATCH --output=run_fit_ldpred2auto_sample_filter_ukb_caucasian_white_british_unrel_by_chr_and_trait.%j.out
#SBATCH --error=run_fit_ldpred2auto_sample_filter_ukb_caucasian_white_british_unrel_by_chr_and_trait.%j.err

###Load module
source /opt/intel/oneapi/mkl/2023.0.0/env/vars.sh intel64
module load R/4.2.3

###Set variables
export MKL_NUM_THREADS=1
export OMP_NUM_THREADS=1

###Get chromosome number from array task id
CHR="${SLURM_ARRAY_TASK_ID}"

###Set some variables
H2_INIT=0.023
DATA_ID=${1}
TRAIT=${2}

###Run R script
echo "Start analyzing trait ${TRAIT} and chromosome ${CHR} ..."

Rscript ../code/fit_ldpred2auto_ukb.R --sumstats ukb_caucasian_white_british_unrel_100000_chr"${CHR}"_sumstats \
                                      --LD_matrix ukb_caucasian_white_british_unrel_100000_chr"${CHR}"_LD \
                                      --n 100000 \
                                      --h2_init "${H2_INIT}" \
                                      --burn_in 500 \
                                      --num_iter 500 \
                                      --trait ${TRAIT} \
                                      --ncores 4 \
                                      --data_id ${DATA_ID} \
                                      --outprefix ukb_caucasian_white_british_unrel_100000_chr"${CHR}"

echo "Finished analyzing trait ${TRAIT} and chromosome ${CHR}."

###Unload module
module unload R/4.2.3

