#!/bin/bash

#SBATCH --job-name=sample_filter_ukb_caucasian_white_british_unrel
#SBATCH --cpus-per-task=36
#SBATCH --partition=fm-bigmem-1,fm-bigmem-2,fm-bigmem-3
#SBATCH --time=4:00:00
#SBATCH --mem=100G
#SBATCH --output=sample_filter_ukb_caucasian_white_british_unrel.%j.out
#SBATCH --error=sample_filter_ukb_caucasian_white_british_unrel.%j.err

/data2/morgante_lab/fabiom/software/plink19 --bfile /data2/morgante_lab/data/ukbiobank/genotypes/array/ukb22418_all_auto_b0_v2_s488243_caucasian_white_british_unrel \
					    --keep-allele-order \
					    --remove /data2/morgante_lab/data/ukbiobank/ind_to_remove/w62347_all_plink.txt \
                                            --maf 0.01 --mac 5 --geno 0.1 --hwe 1e-10 \
                                            --thin-indiv-count 105000 \
                                            --make-bed \
                                            --out ../data/genotypes/ukb_caucasian_white_british_unrel_100000 \
					    --seed 123 \
                                            --threads 36 

