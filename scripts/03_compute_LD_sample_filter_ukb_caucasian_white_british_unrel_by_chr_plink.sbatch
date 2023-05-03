#!/bin/bash

#SBATCH --job-name=compute_LD_sample_filter_ukb_caucasian_white_british_unrel_by_chr_plink
#SBATCH --cpus-per-task=16
#SBATCH --partition=fm-bigmem-1,fm-bigmem-2,fm-bigmem-3
#SBATCH --time=05:00:00
#SBATCH --mem=50G
#SBATCH --output=compute_LD_sample_filter_ukb_caucasian_white_british_unrel_by_chr_plink.%j.out
#SBATCH --error=compute_LD_sample_filter_ukb_caucasian_white_british_unrel_by_chr_plink.%j.err

###Get chromosome number from array task id
CHR="${SLURM_ARRAY_TASK_ID}"

###Run plink
echo "Start analyzing chromosome ${CHR} ..."

/data2/morgante_lab/fabiom/software/plink19 --bfile /data2/morgante_lab/data/ukbiobank/genotypes/array/ukb22418_all_auto_b0_v2_s488243_caucasian_white_british_unrel \
                                            --keep-allele-order \
                                            --remove ../data/LD_matrices/inds_to_remove.txt \
                                            --extract ../data/genotypes/ukb_caucasian_white_british_unrel_100000.snplist \
                                            --chr ${CHR} \
                                            --r bin 'yes-really' \
                                            --out ../data/LD_matrices/ukb_caucasian_white_british_unrel_chr"${CHR}"_LD \
                                            --threads 16

echo "Finished analyzing chromosome ${CHR}."
