#!/bin/bash

#SBATCH --job-name=compute_LD_scores_1000G_for_MTAG
#SBATCH --cpus-per-task=4
#SBATCH --partition=fm-bigmem-1,fm-bigmem-2,fm-bigmem-3,fm-bigmem-4
#SBATCH --time=24:00:00
#SBATCH --mem=10G
#SBATCH --output=compute_LD_scores_1000G_for_MTAG.%j.out
#SBATCH --error=compute_LD_scores_1000G_for_MTAG.%j.err

###Get chromosome number from array task id
CHR="${SLURM_ARRAY_TASK_ID}"

source /data2/morgante_lab/fabiom/software/miniconda3/etc/profile.d/conda.sh
conda activate mtag

python /data2/morgante_lab/fabiom/software/ldsc/ldsc.py \
	--bfile "../data/genotypes/1000G_phase3_eur_for_LD_genetic_dist_chr${CHR}" \
	--l2 \
	--ld-wind-cm 1 \
	--out ../data/LD_scores/simulations_1000G/${CHR}

conda deactivate

