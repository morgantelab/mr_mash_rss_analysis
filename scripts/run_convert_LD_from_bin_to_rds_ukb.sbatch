#!/bin/bash

#SBATCH --job-name=run_convert_LD_from_bin_to_rds_ukb
#SBATCH --cpus-per-task=1
#SBATCH --partition=fm-bigmem-1,fm-bigmem-2,fm-bigmem-3,fm-bigmem-4
#SBATCH --time=06:00:00
#SBATCH --mem=175G
#SBATCH --output=run_convert_LD_from_bin_to_rds_ukb.%j.out
#SBATCH --error=run_convert_LD_from_bin_to_rds_ukb.%j.err

###Load module
source /opt/intel/oneapi/mkl/2023.0.0/env/vars.sh intel64
module load R/4.2.3

###Set variables
export MKL_NUM_THREADS=1

###Get chromosome number from array task id
CHR="${SLURM_ARRAY_TASK_ID}"

Rscript ../code/convert_LD_from_bin_to_rds.R --sumstats ../output/summary_statistics/ukb_caucasian_white_british_unrel_100000_equal_effects_indep_resid_chr"${CHR}"_sumstats_1.rds \
                                             --LD_matrix ../data/LD_matrices/ukb_caucasian_white_british_unrel_chr"${CHR}"_LD_sparse.ld.bin \
                                             --output ../data/LD_matrices/ukb_caucasian_white_british_unrel_chr"${CHR}"_LD_sparse.rds

