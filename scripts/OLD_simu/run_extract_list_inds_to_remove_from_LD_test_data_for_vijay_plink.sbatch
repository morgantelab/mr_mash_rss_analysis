#!/bin/bash

#SBATCH --job-name=run_extract_list_inds_to_remove_from_LD_test_data_for_vijay_plink
#SBATCH --cpus-per-task=1
#SBATCH --partition=fm-bigmem-1,fm-bigmem-2,fm-bigmem-3
#SBATCH --time=00:30:00
#SBATCH --mem=2G
#SBATCH --output=run_extract_list_inds_to_remove_from_LD_test_data_for_vijay_plink.%j.out
#SBATCH --error=run_extract_list_inds_to_remove_from_LD_test_data_for_vijay_plink.%j.err

###Load module
source /opt/intel/oneapi/mkl/2023.0.0/env/vars.sh intel64
module load R/4.2.3

###Set variables
export MKL_NUM_THREADS=1
export OMP_NUM_THREADS=1

###Run script
Rscript ../code/extract_list_inds_to_remove_from_LD_ukb.R --geno ../data/genotypes/test_data_for_vijay.fam \
                                                          --withdrawn_inds /data2/morgante_lab/data/ukbiobank/ind_to_remove/w62347_all_plink.txt \
                                                          --output ../data/LD_matrices/inds_to_remove_test_data_for_vijay.txt

###Unload module
module unload R/4.2.3
