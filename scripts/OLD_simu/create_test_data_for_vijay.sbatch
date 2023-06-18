#!/bin/bash

#SBATCH --job-name=create_test_data_for_vijay
#SBATCH --cpus-per-task=1
#SBATCH --partition=fm-bigmem-1,fm-bigmem-2,fm-bigmem-3
#SBATCH --time=4:00:00
#SBATCH --mem=50G
#SBATCH --output=create_test_data_for_vijay.%j.out
#SBATCH --error=create_test_data_for_vijay.%j.err

/data2/morgante_lab/fabiom/software/plink19 --bfile /data2/morgante_lab/ukbiobank_projects/mr_mash_rss/data/genotypes/ukb_caucasian_white_british_unrel_100000 \
                                            --keep-allele-order \
                                            --thin-indiv-count 1000 \
                                            --thin-count 5000 \
                                            --write-snplist \
                                            --make-bed \
                                            --out /data2/morgante_lab/ukbiobank_projects/mr_mash_rss/data/genotypes/test_data_for_vijay \
                                            --seed 1 \
                                            --threads 1 

