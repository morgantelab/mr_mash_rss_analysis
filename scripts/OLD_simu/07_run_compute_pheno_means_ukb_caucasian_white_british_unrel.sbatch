#!/bin/bash

#SBATCH --job-name=run_compute_pheno_means_ukb_caucasian_white_british_unrel_by_chr
#SBATCH --cpus-per-task=1
#SBATCH --partition=fm-bigmem-1,fm-bigmem-2,fm-bigmem-3
#SBATCH --time=00:30:00
#SBATCH --mem=2G
#SBATCH --output=run_compute_pheno_means_ukb_caucasian_white_british_unrel.%j.out
#SBATCH --error=run_compute_pheno_means_ukb_caucasian_white_british_unrel.%j.err

###Set some variable
DATA_ID=${1}

###Load module
source /opt/intel/oneapi/mkl/2023.0.0/env/vars.sh intel64
module load R/4.2.3

###Set variables
export MKL_NUM_THREADS=1
export OMP_NUM_THREADS=1

###Run R script
Rscript ../code/compute_phenotype_means_ukb.R --pheno ../data/phenotypes/simulated/ukb_caucasian_white_british_unrel_100000_pheno_"${DATA_ID}".rds \
                                             --test_ids ../data/phenotypes/simulated/ukb_caucasian_white_british_unrel_100000_test_ids_"${DATA_ID}".txt \
                                             --output ../output/misc/ukb_caucasian_white_british_unrel_100000_pheno_means_"${DATA_ID}".rds \
                                             --seed ${DATA_ID}

###Unload module
module unload R/4.2.3

