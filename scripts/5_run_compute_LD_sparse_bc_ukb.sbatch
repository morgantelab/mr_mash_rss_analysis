#!/bin/bash

#SBATCH --job-name=run_compute_LD_sparse_bc_by_fold_ukb
#SBATCH --cpus-per-task=16
#SBATCH --partition=fm-bigmem-1,fm-bigmem-2,fm-bigmem-3
#SBATCH --time=48:00:00
#SBATCH --mem=150G
#SBATCH --output=run_compute_LD_sparse_bc_by_fold_ukb.%j.out
#SBATCH --error=run_compute_LD_sparse_bc_by_fold_ukb.%j.err

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

Rscript ../code/compute_LD_matrix_bigsnp_obj_ukb.R --geno /scratch1/fabiom/ukb_bc_geno_imp_HM3.rds \
                                                   --chr ${CHR} \
                                                   --samples_to_remove ../data/misc/ukb_cleaned_bc_ind_ids_plink_"${FOLD}".txt \
                                                   --sparse TRUE \
                                                   --genetic_map ../data/misc/1000-genomes-genetic-maps/interpolated_OMNI \
                                                   --window_size 0.003 \
                                                   --impute_missing FALSE \
                                                   --seed 1 \
                                                   --output ../data/LD_matrices/ukb_bc_chr"${CHR}"_LD_sparse_"${FOLD}".ld.bin \
                                                   --ncores 8

echo "Finished analyzing chromosome ${CHR} and fold ${FOLD}."

###Unload module
module unload R/4.2.3

