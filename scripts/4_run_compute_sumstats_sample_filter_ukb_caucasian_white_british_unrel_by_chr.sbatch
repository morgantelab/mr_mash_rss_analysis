#!/bin/bash

#SBATCH --job-name=run_compute_sumstats_sample_filter_ukb_caucasian_white_british_unrel_by_chr
#SBATCH --cpus-per-task=1
#SBATCH --partition=fm-bigmem-1,fm-bigmem-2,fm-bigmem-3
#SBATCH --time=24:00:00
#SBATCH --mem=100G
#SBATCH --output=run_compute_sumstats_sample_filter_ukb_caucasian_white_british_unrel_by_chr.%j.out
#SBATCH --error=run_compute_sumstats_sample_filter_ukb_caucasian_white_british_unrel_by_chr.%j.err

###Load module
module load gsl/2.7
source /opt/intel/oneapi/mkl/2023.0.0/env/vars.sh intel64
module load R/4.2.3

###Set variables
export MKL_NUM_THREADS=1
export OMP_NUM_THREADS=1

###Get chromosome number from array task id
CHR="${SLURM_ARRAY_TASK_ID}"

###Run R script
echo "Start analyzing chromosome ${CHR} ..."

Rscript ../code/compute_sumstats_ukb.R --pheno simulated/ukb_caucasian_white_british_unrel_100000_pheno \
                                       --geno ukb_caucasian_white_british_unrel_100000 \
                                       --test_ids simulated/ukb_caucasian_white_british_unrel_100000_test_ids \
                                       --outprefix ukb_caucasian_white_british_unrel_100000 \
                                       --data_id 1 \
                                       --chr ${CHR} \
                                       --standardize FALSE \
                                       --impute_missing TRUE \
                                       --ncores 1

echo "Finished analyzing chromosome ${CHR}."

###Unload module
module unload R/4.2.3
module unload gsl/2.7

