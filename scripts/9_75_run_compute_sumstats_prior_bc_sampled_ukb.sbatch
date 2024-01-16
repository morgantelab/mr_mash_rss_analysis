#!/bin/bash

#SBATCH --job-name=run_compute_sumstats_prior_bc_sampled_ukb
#SBATCH --cpus-per-task=16
#SBATCH --partition=fm-bigmem-1,fm-bigmem-2,fm-bigmem-3
#SBATCH --time=12:00:00
#SBATCH --mem=100G
#SBATCH --output=run_compute_sumstats_prior_bc_sampled_ukb.%j.out
#SBATCH --error=run_compute_sumstats_prior_bc_sampled_ukb.%j.err

###Load module
source /opt/intel/oneapi/mkl/2023.0.0/env/vars.sh intel64
module load R/4.2.3

###Get chromosome number
CHR=0

###Get fold from arguments
FOLD=${1}

###Run R script
echo "Start analyzing chromosome ${CHR} and fold ${FOLD} ..."

Rscript ../code/compute_sumstats_bc_sampled_prior_ukb.R --fold ${FOLD} \
                                                        --chr ${CHR} \
                                                        --maf 0.01 \
                                                        --standardize FALSE \
                                                        --impute_missing FALSE \
                                                        --normalize TRUE \
                                                        --ncores 16

echo "Finished analyzing chromosome ${CHR} and fold ${FOLD}."

###Unload module
module unload R/4.2.3

