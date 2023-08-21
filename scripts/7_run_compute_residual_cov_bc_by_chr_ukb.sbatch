#!/bin/bash

#SBATCH --job-name=run_compute_residual_cov_bc_by_chr_ukb
#SBATCH --cpus-per-task=1
#SBATCH --partition=fm-bigmem-1,fm-bigmem-2,fm-bigmem-3
#SBATCH --time=01:00:00
#SBATCH --mem=4G
#SBATCH --output=run_compute_residual_cov_bc_by_chr_ukb.%j.out
#SBATCH --error=run_compute_residual_cov_bc_by_chr_ukb.%j.err

###Get chromosome number from array task id
CHR="${SLURM_ARRAY_TASK_ID}"

###Set some variable
FOLD=${1}

###Load module
module load gsl/2.7
source /opt/intel/oneapi/mkl/2023.0.0/env/vars.sh intel64
module load R/4.2.3

###Run R script
echo "Start analyzing chromosome ${CHR} and fold ${FOLD} ..."

Rscript ../code/compute_residual_cov_all_chr_ukb.R --sumstats_prefix ../output/summary_statistics/ukb_bc_chr \
                                                   --sumstats_suffix _sumstats_"${FOLD}".rds \
                                                   --chr ${CHR} \
                                                   --output ../output/misc/ukb_bc_chr"${CHR}"_residual_cov_"${FOLD}".rds \
                                                   --seed ${FOLD}

echo "Finished analyzing chromosome ${CHR} and fold ${FOLD}."

###Unload module
module unload R/4.2.3
module load gsl/2.7

