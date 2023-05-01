#!/bin/bash

#SBATCH --job-name=run_split_train_test_sample_filter_ukb_caucasian_white_british_unrel_list_plink
#SBATCH --cpus-per-task=1
#SBATCH --partition=fm-bigmem-1,fm-bigmem-2,fm-bigmem-3
#SBATCH --time=00:10:00
#SBATCH --mem=2G
#SBATCH --output=run_split_train_test_sample_filter_ukb_caucasian_white_british_unrel_list_plink.%j.out
#SBATCH --error=run_split_train_test_sample_filter_ukb_caucasian_white_british_unrel_list_plink.%j.err

###Set some variable
DATA_ID=${1}

###Load module
source /opt/intel/oneapi/mkl/2023.0.0/env/vars.sh intel64
module load R/4.2.3

###Run script
Rscript ../code/split_train_test_ukb_list_plink.R --input ukb_caucasian_white_british_unrel_100000 \
						  --n_train 100000 \
                                                  --data_id ${DATA_ID}

###Unload module
module unload R/4.2.3

