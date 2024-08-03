#!/bin/bash

#SBATCH --job-name=run_compute_LD_block_sparse_1000G_for_sims_by_chr_bigsnpr
#SBATCH --cpus-per-task=8
#SBATCH --partition=fm-bigmem-1,fm-bigmem-2,fm-bigmem-3,fm-bigmem-4
#SBATCH --time=02:00:00
#SBATCH --mem=100G
#SBATCH --output=run_compute_LD_block_sparse_1000G_for_sims_by_chr_bigsnpr.%j.out
#SBATCH --error=run_compute_LD_block_sparse_1000G_for_sims_by_chr_bigsnpr.%j.err

###Load module
source /opt/intel/oneapi/mkl/2023.0.0/env/vars.sh intel64
module load R/4.2.3

###Set variables
export MKL_NUM_THREADS=8
export OMP_NUM_THREADS=8

###Get chromosome number from array task id
CHR="${SLURM_ARRAY_TASK_ID}"

###Run R script
echo "Start analyzing chromosome ${CHR} ..."

Rscript ../code/compute_LD_matrix_blocks_1000G.R --geno ../data/genotypes/1000G_phase3_eur_for_LD_genetic_dist_chr"${CHR}".bed \
                                                 --thr_maf 0.05 \
                                                 --thr_r2_eigen 0.8 \
                                                 --blocks /data2/morgante_lab/data/1000G/ldetect-data/EUR/fourier_ls-chr"${CHR}".bed \
                                                 --impute_missing FALSE \
                                                 --seed 1 \
                                                 --temp_dir /scratch3/fabiom/tmp \
                                                 --output ../data/LD_matrices/1000G_phase3_eur_chr"${CHR}"_LD_block_sparse \
                                                 --ncores 1

echo "Finished analyzing chromosome ${CHR}."

###Unload module
module unload R/4.2.3

