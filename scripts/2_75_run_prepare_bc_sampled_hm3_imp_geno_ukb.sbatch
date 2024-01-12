#!/bin/bash

#SBATCH --job-name=run_prepare_bc_sampled_hm3_imp_geno_ukb
#SBATCH --cpus-per-task=22
#SBATCH --partition=fm-bigmem-1,fm-bigmem-2,fm-bigmem-3,bigmem
#SBATCH --time=12:00:00
#SBATCH --mem=250G
#SBATCH --output=run_prepare_bc_sampled_hm3_imp_geno_ukb.%j.out
#SBATCH --error=run_prepare_bc_sampled_hm3_imp_geno_ukb.%j.err

###Load module
source /opt/intel/oneapi/mkl/2023.0.0/env/vars.sh intel64
module load R/4.2.3

###Run script
Rscript ../code/prepare_bc_sampled_hm3_imp_geno_ukb.R

###Unload module
module unload R/4.2.3
