#!/bin/bash

#SBATCH --job-name=run_compute_phenotypic_cov_ukb_caucasian_white_british_unrel
#SBATCH --cpus-per-task=8
#SBATCH --partition=fm-bigmem-1,fm-bigmem-2,fm-bigmem-3
#SBATCH --time=02:00:00
#SBATCH --mem=6G
#SBATCH --output=run_compute_phenotypic_cov_ukb_caucasian_white_british_unrel.%j.out
#SBATCH --error=run_compute_phenotypic_cov_ukb_caucasian_white_british_unrel.%j.err

###Load module
module load gsl/2.7
source /opt/intel/oneapi/mkl/2023.0.0/env/vars.sh intel64
module load R/4.2.3

###Run R script
Rscript ../code/compute_phenotypic_cov_ukb_sim.R --input ukb_caucasian_white_british_unrel_100000 \
                                                 --test_ids simulated/ukb_caucasian_white_british_unrel_100000_test_ids \
                                                 --data_id 1

###Unload module
module unload R/4.2.3
module unload gsl/2.7
