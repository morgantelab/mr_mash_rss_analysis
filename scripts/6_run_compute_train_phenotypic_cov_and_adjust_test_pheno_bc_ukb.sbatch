#!/bin/bash

#SBATCH --job-name=run_compute_train_phenotypic_cov_and_adjust_test_pheno_bc_ukb
#SBATCH --cpus-per-task=1
#SBATCH --partition=fm-bigmem-1,fm-bigmem-2,fm-bigmem-3
#SBATCH --time=02:00:00
#SBATCH --mem=15G
#SBATCH --output=run_compute_train_phenotypic_cov_and_adjust_test_pheno_bc_ukb.%j.out
#SBATCH --error=run_compute_train_phenotypic_cov_and_adjust_test_pheno_bc_ukb.%j.err

###Load module
source /opt/intel/oneapi/mkl/2023.0.0/env/vars.sh intel64
module load R/4.2.3

###Run script
Rscript ../code/compute_train_phenotypic_cov_and_adjust_test_pheno_BP_ukb.R --normalize FALSE

###Unload module
module unload R/4.2.3
