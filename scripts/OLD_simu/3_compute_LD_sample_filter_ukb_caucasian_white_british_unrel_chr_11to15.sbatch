#!/bin/bash

#SBATCH --job-name=compute_LD_sample_filter_ukb_caucasian_white_british_unrel_chr_11to15
#SBATCH --cpus-per-task=16
#SBATCH --partition=fm-bigmem-1,fm-bigmem-2,fm-bigmem-3
#SBATCH --time=12:00:00
#SBATCH --mem=30G
#SBATCH --output=compute_LD_sample_filter_ukb_caucasian_white_british_unrel_chr_11to15.%j.out
#SBATCH --error=compute_LD_sample_filter_ukb_caucasian_white_british_unrel_chr_11to15.%j.err

data_id=1

for i in {11..15}
do
   /data2/morgante_lab/fabiom/software/plink19 --bfile ../data/genotypes/ukb_caucasian_white_british_unrel_100000 \
                                               --keep-allele-order \
                                               --remove ../data/phenotypes/simulated/ukb_caucasian_white_british_unrel_100000_test_ids_"${data_id}".txt \
                                               --chr ${i} \
                                               --r bin 'yes-really' \
                                               --out ../data/LD_matrices/ukb_caucasian_white_british_unrel_100000_chr"${i}"_LD_"${data_id}" \
                                               --threads 16
done

