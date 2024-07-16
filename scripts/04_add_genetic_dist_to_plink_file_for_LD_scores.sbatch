#!/bin/bash

#SBATCH --job-name=add_genetic_dist_to_plink_file_for_LD_scores
#SBATCH --cpus-per-task=1
#SBATCH --partition=fm-bigmem-1,fm-bigmem-2,fm-bigmem-3
#SBATCH --time=01:00:00
#SBATCH --mem=10G
#SBATCH --output=add_genetic_dist_to_plink_file.%j.out
#SBATCH --error=add_genetic_dist_to_plink_file.%j.err

# wget https://www.dropbox.com/s/slchsd0uyd4hii8/genetic_map_b37.zip
# unzip genetic_map_b37.zip

###Annotate genetic distance
/data2/morgante_lab/fabiom/software/plink19 --bfile /data2/morgante_lab/data/ukbiobank/genotypes/array/ukb22418_all_auto_b0_v2_s488243_caucasian_white_british_unrel \
                                            --keep-allele-order \
                                            --remove ../data/LD_matrices/inds_to_remove.txt \
                                            --extract ../data/genotypes/ukb_caucasian_white_british_unrel_100000.snplist \
                                            --cm-map ../data/misc/1000-genomes-genetic-map-ldak/genetic_map_b37/genetic_map_chr@_combined_b37.txt \
                                            --make-bed \
                                            --out ../data/genotypes/ukb_caucasian_white_british_unrel_for_LD_genetic_dist



###Split data into chromosomes
for CHR in {1..22};
do
  /data2/morgante_lab/fabiom/software/plink19 --bfile ../data/genotypes/ukb_caucasian_white_british_unrel_for_LD_genetic_dist \
                                              --keep-allele-order \
                                              --chr ${CHR} \
                                              --make-bed \
                                              --out "../data/genotypes/ukb_caucasian_white_british_unrel_for_LD_genetic_dist_chr${CHR}"
done

rm ../data/genotypes/ukb_caucasian_white_british_unrel_for_LD_genetic_dist.bed ../data/genotypes/ukb_caucasian_white_british_unrel_for_LD_genetic_dist.bim ../data/genotypes/ukb_caucasian_white_british_unrel_for_LD_genetic_dist.fam
