#!/bin/bash

#SBATCH --job-name=run_compute_LD_sparse_bc_ukb
#SBATCH --cpus-per-task=8
#SBATCH --partition=fm-bigmem-1,fm-bigmem-2,fm-bigmem-3
#SBATCH --time=48:00:00
#SBATCH --mem=100G
#SBATCH --output=run_compute_LD_sparse_bc_ukb.%j.out
#SBATCH --error=run_compute_LD_sparse_bc_ukb.%j.err

###Load module
source /opt/intel/oneapi/mkl/2023.0.0/env/vars.sh intel64
module load R/4.2.3

###Set variables
export MKL_NUM_THREADS=1

###Get chromosome number from array task id
CHR="${SLURM_ARRAY_TASK_ID}"

###Run R script
echo "Start analyzing chromosome ${CHR} ..."

Rscript ../code/compute_LD_matrix_bigsnp_obj_ukb.R --geno /scratch1/fabiom/ukb_bc_geno_imp_HM3.rds \
                                                   --chr ${CHR} \
                                                   --sparse TRUE \
                                                   --genetic_map ../data/misc/1000-genomes-genetic-maps/interpolated_OMNI \
                                                   --window_size 0.003 \
                                                   --impute_missing FALSE \
                                                   --seed 1 \
                                                   --output ../data/LD_matrices/ukb_bc_chr"${CHR}"_LD_sparse.ld.bin \
                                                   --ncores 8

echo "Finished analyzing chromosome ${CHR}."

###Unload module
module unload R/4.2.3

