#!/bin/bash

#SBATCH --job-name=run_compute_geno_means_ukb_caucasian_white_british_unrel_by_chr
#SBATCH --cpus-per-task=1
#SBATCH --partition=fm-bigmem-1,fm-bigmem-2,fm-bigmem-3
#SBATCH --time=24:00:00
#SBATCH --mem=100G
#SBATCH --output=run_compute_geno_means_ukb_caucasian_white_british_unrel_by_chr.%j.out
#SBATCH --error=run_compute_geno_means_ukb_caucasian_white_british_unrel_by_chr.%j.err

###Set some variable
DATA_ID=${1}

###Load module
source /opt/intel/oneapi/mkl/2023.0.0/env/vars.sh intel64
module load R/4.2.3

###Set variables
export MKL_NUM_THREADS=1
export OMP_NUM_THREADS=1

###Get chromosome number from array task id
CHR="${SLURM_ARRAY_TASK_ID}"

###Run R script
echo "Start analyzing chromosome ${CHR} ..."

Rscript ../code/compute_genotype_means_ukb.R --geno ../data/genotypes/ukb_caucasian_white_british_unrel_100000.bed \
                                             --test_ids ../data/phenotypes/simulated/ukb_caucasian_white_british_unrel_100000_test_ids_"${DATA_ID}".txt \
                                             --output ../output/misc/ukb_caucasian_white_british_unrel_100000_chr"${CHR}"_geno_means_"${DATA_ID}".rds \
                                             --seed ${DATA_ID} \
                                             --chr ${CHR} \
                                             --impute_missing TRUE \
                                             --temp_dir /data2/morgante_lab/fabiom/tmp \
                                             --ncores 1

echo "Finished analyzing chromosome ${CHR}."

###Unload module
module unload R/4.2.3

