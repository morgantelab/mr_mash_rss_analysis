#!/bin/bash

#SBATCH --job-name=run_compute_sumstats_BP_by_chr_and_fold
#SBATCH --cpus-per-task=1
#SBATCH --partition=fm-bigmem-1,fm-bigmem-2,fm-bigmem-3
#SBATCH --time=12:00:00
#SBATCH --mem=20G
#SBATCH --output=run_compute_sumstats_BP_by_chr_and_fold.%j.out
#SBATCH --error=run_compute_sumstats_BP_by_chr_and_fold.%j.err

###Load module
source /opt/intel/oneapi/mkl/2023.0.0/env/vars.sh intel64
module load R/4.2.3

###Get chromosome number from array task id
CHR="${SLURM_ARRAY_TASK_ID}"

###Get fold from arguments
FOLD=${1}

###Run R script
echo "Start analyzing chromosome ${CHR} and fold ${FOLD} ..."

Rscript ../code/compute_sumstats_BP_ukb.R --fold ${FOLD} \
                                          --chr ${CHR} \
                                          --standardize FALSE \
                                          --impute_missing FALSE \
                                          --normalize TRUE \
                                          --ncores 3

echo "Finished analyzing chromosome ${CHR} and fold ${FOLD}."

###Unload module
module unload R/4.2.3

