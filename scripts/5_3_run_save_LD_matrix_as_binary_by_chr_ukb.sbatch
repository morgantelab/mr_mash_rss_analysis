#!/bin/bash

#SBATCH --job-name=run_save_LD_matrix_as_binary_by_chr_ukb
#SBATCH --cpus-per-task=10
#SBATCH --partition=fm-bigmem-1,fm-bigmem-2,fm-bigmem-3
#SBATCH --time=01:00:00
#SBATCH --mem=150G
#SBATCH --output=run_save_LD_matrix_as_binary_by_chr_ukb.%j.out
#SBATCH --error=run_save_LD_matrix_as_binary_by_chr_ukb.%j.err

###Load module
source /opt/intel/oneapi/mkl/2023.0.0/env/vars.sh intel64
module load R/4.2.3

###Set variables
export MKL_NUM_THREADS=1

###Get chromosome number from array task id
CHR="${SLURM_ARRAY_TASK_ID}"

###Get fold from arguments
FOLD=${1}

###Run R script
echo "Start analyzing chromosome ${CHR} and fold ${FOLD} ..."

Rscript ../code/save_LD_matrix_as_binary_ukb.R --input /data2/morgante_lab/ukbiobank_projects/mr_mash_rss/tmp/ukb_bc_chr"${CHR}"_LD_"${FOLD}".ld \
                                               --output ../data/LD_matrices/ukb_bc_chr"${CHR}"_LD_"${FOLD}".ld.bin \
                                               --ncores 10

echo "Finished analyzing chromosome ${CHR} and fold ${FOLD}."

###Unload module
module unload R/4.2.3
