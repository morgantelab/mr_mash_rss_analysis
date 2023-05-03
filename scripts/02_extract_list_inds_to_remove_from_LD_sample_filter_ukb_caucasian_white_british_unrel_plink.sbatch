#!/bin/bash

#SBATCH --job-name=extract_list_inds_to_remove_from_LD_sample_filter_ukb_caucasian_white_british_unrel_plink
#SBATCH --cpus-per-task=1
#SBATCH --partition=fm-bigmem-1,fm-bigmem-2,fm-bigmem-3
#SBATCH --time=00:30:00
#SBATCH --mem=1G
#SBATCH --output=extract_list_inds_to_remove_from_LD_sample_filter_ukb_caucasian_white_british_unrel_plink.%j.out
#SBATCH --error=extract_list_inds_to_remove_from_LD_sample_filter_ukb_caucasian_white_british_unrel_plink.%j.err

###Concatenate list of individuals to remove
cut -d" " -f1-2 ../data/genotypes/ukb_caucasian_white_british_unrel_100000.fam | cat - /data2/morgante_lab/data/ukbiobank/ind_to_remove/w62347_all_plink.txt > ../data/LD_matrices/inds_to_remove.txt

