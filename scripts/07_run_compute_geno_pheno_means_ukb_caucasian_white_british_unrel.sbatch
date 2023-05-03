#!/bin/bash

#SBATCH --job-name=run_compute_geno_pheno_means_ukb_caucasian_white_british_unrel_by_chr
#SBATCH --cpus-per-task=16
#SBATCH --partition=fm-bigmem-1,fm-bigmem-2,fm-bigmem-3
#SBATCH --time=24:00:00
#SBATCH --mem=100G
#SBATCH --output=run_compute_geno_pheno_means_ukb_caucasian_white_british_unrel_by_chr.%j.out
#SBATCH --error=run_compute_geno_pheno_means_ukb_caucasian_white_british_unrel_by_chr.%j.err

###Set some variable
DATA_ID=${1}

###Load module
source /opt/intel/oneapi/mkl/2023.0.0/env/vars.sh intel64
module load R/4.2.3

###Set variables
export MKL_NUM_THREADS=1

###Run R script
Rscript ../code/compute_genotype_phenotype_means_ukb.R --pheno simulated/ukb_caucasian_white_british_unrel_100000_pheno \
                                                       --geno ukb_caucasian_white_british_unrel_100000 \
                                                       --test_ids simulated/ukb_caucasian_white_british_unrel_100000_test_ids \
                                                       --outprefix ukb_caucasian_white_british_unrel_100000 \
                                                       --data_id ${DATA_ID} \
                                                       --chr 1:22 \
                                                       --impute_missing TRUE \
                                                       --ncores 16

###Unload module
module unload R/4.2.3

