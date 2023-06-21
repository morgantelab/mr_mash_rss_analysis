#!/bin/bash

#SBATCH --job-name=run_extract_lists_for_LD_BP_ukb
#SBATCH --cpus-per-task=1
#SBATCH --partition=fm-bigmem-1,fm-bigmem-2,fm-bigmem-3
#SBATCH --time=00:30:00
#SBATCH --mem=4G
#SBATCH --output=run_extract_lists_for_LD_BP_ukb.%j.out
#SBATCH --error=run_extract_lists_for_LD_BP_ukb.%j.err

###Load module
source /opt/intel/oneapi/mkl/2023.0.0/env/vars.sh intel64
module load R/4.2.3

###Run script
Rscript ../code/extract_lists_for_LD_BP_ukb.R

###Unload module
module unload R/4.2.3
