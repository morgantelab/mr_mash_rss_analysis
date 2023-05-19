#!/bin/bash

#SBATCH --job-name=run_extract_list_inds_to_remove_from_LD_sample_filter_ukb_caucasian_white_british_unrel_plink
#SBATCH --cpus-per-task=1
#SBATCH --partition=fm-bigmem-1,fm-bigmem-2,fm-bigmem-3
#SBATCH --time=00:30:00
#SBATCH --mem=2G
#SBATCH --output=run_extract_list_inds_to_remove_from_LD_sample_filter_ukb_caucasian_white_british_unrel_plink.%j.out
#SBATCH --error=run_extract_list_inds_to_remove_from_LD_sample_filter_ukb_caucasian_white_british_unrel_plink.%j.err

###Load module
source /opt/intel/oneapi/mkl/2023.0.0/env/vars.sh intel64
module load R/4.2.3

###Set variables
export MKL_NUM_THREADS=1
export OMP_NUM_THREADS=1

###Run script
Rscript ../code/extract_list_inds_to_remove_from_LD_ukb.R --geno ../data/genotypes/ukb_caucasian_white_british_unrel_100000.fam \
				                          --withdrawn_inds /data2/morgante_lab/data/ukbiobank/ind_to_remove/w62347_all_plink.txt \
                                --output ../data/LD_matrices/inds_to_remove.txt

###Unload module
module unload R/4.2.3
