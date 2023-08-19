#!/bin/bash

#SBATCH --job-name=run_prepare_tiezzi_BP_pheno_ukb
#SBATCH --cpus-per-task=1
#SBATCH --partition=fm-bigmem-1,fm-bigmem-2,fm-bigmem-3
#SBATCH --time=01:00:00
#SBATCH --mem=100G
#SBATCH --output=run_prepare_tiezzi_BP_pheno_ukb.%j.out
#SBATCH --error=run_prepare_tiezzi_BP_pheno_ukb.%j.err

###Load module
source /opt/intel/oneapi/mkl/2023.0.0/env/vars.sh intel64
module load R/4.2.3

###Set variables
export MKL_NUM_THREADS=1
export OMP_NUM_THREADS=1

###Run script
Rscript ../code/prepare_tiezzi_BP_pheno_ukb.R

###Unload module
module unload R/4.2.3
