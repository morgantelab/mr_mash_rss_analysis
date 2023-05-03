#!/bin/bash

#SBATCH --job-name=run_sim_pheno_sample_filter_ukb_caucasian_white_british_unrel
#SBATCH --cpus-per-task=16
#SBATCH --partition=fm-bigmem-1,fm-bigmem-2,fm-bigmem-3
#SBATCH --time=04:00:00
#SBATCH --mem=200G
#SBATCH --output=run_sim_pheno_sample_filter_ukb_caucasian_white_british_unrel.%j.out
#SBATCH --error=run_sim_pheno_sample_filter_ukb_caucasian_white_british_unrel.%j.err


###Set some variable
DATA_ID=${1}

###Load module
module load gsl/2.7
source /opt/intel/oneapi/mkl/2023.0.0/env/vars.sh intel64
module load R/4.2.3

###Set variables
export MKL_NUM_THREADS=16
export OMP_NUM_THREADS=16

###Run script
Rscript ../code/sim_pheno_ukb.R --geno ukb_caucasian_white_british_unrel_100000 \
				--p_causal 2000 --r 5 --pve 0.3 --B_cor 1 \
				--V_cor 0 --data_id ${DATA_ID}

###Unload module
module unload R/4.2.3
module unload gsl/2.7

