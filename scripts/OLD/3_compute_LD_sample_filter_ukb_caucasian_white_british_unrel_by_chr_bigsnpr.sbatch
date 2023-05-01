#!/bin/bash

#SBATCH --job-name=compute_LD_sample_filter_ukb_caucasian_white_british_unrel_by_chr
#SBATCH --cpus-per-task=4
#SBATCH --partition=fm-bigmem-1,fm-bigmem-2,fm-bigmem-3
#SBATCH --time=24:00:00
#SBATCH --mem=80G
#SBATCH --output=compute_LD_sample_filter_ukb_caucasian_white_british_unrel_by_chr.%j.out
#SBATCH --error=compute_LD_sample_filter_ukb_caucasian_white_british_unrel_by_chr.%j.err

###Load module
source /opt/intel/oneapi/mkl/2023.0.0/env/vars.sh intel64
module load R/4.2.3

###Set variables
export MKL_NUM_THREADS=1

###Get chromosome number from array task id
CHR="${SLURM_ARRAY_TASK_ID}"

###Get data_id from arguments
DATA_ID=${1}

###Run R script
echo "Start analyzing chromosome ${CHR} ..."

Rscript ../code/compute_LD_matrix_ukb.R --geno ukb_caucasian_white_british_unrel_100000 \
                                        --test_ids simulated/ukb_caucasian_white_british_unrel_100000_test_ids \
                                        --data_id ${DATA_ID} \
                                        --chr ${CHR} \
                                        --impute_missing TRUE \
                                        --ncores 4

echo "Finished analyzing chromosome ${CHR}."

###Unload module
module unload R/4.2.3

