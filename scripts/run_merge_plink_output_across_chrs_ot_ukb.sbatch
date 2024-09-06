#!/bin/bash

#SBATCH --job-name=run_merge_plink_output_across_chrs_ot_ukb
#SBATCH --cpus-per-task=8
#SBATCH --partition=fm-bigmem-1,fm-bigmem-2,fm-bigmem-3,fm-bigmem-4
#SBATCH --time=02:00:00
#SBATCH --mem=10G
#SBATCH --output=run_merge_plink_output_across_chrs_ot_ukb.%j.out
#SBATCH --error=run_merge_plink_output_across_chrs_ot_ukb.%j.err

###Load module
source /opt/intel/oneapi/mkl/2023.0.0/env/vars.sh intel64
module load R/4.2.3

###Set variables
export MKL_NUM_THREADS=1

for TRAIT in weight waist hip BMI TFM BFP DPa SPa;
do
    Rscript ../code/merge_plink_output_across_chrs_ot_ukb.R --prefix ../output/GWAS_for_regions/ukb_ot_chr \
                                                            --suffix _biallelic_snps_only_recoded_ids_info_06_gwas_for_regions.rsid."${TRAIT}".glm.linear \
       		                                                  --chr 1:22 \
           		                                              --ncores 8 \
               		                                          --output ../output/GWAS_for_regions/ukb_ot_chrAll_biallelic_snps_only_recoded_ids_info_06_gwas_for_regions.rsid."${TRAIT}".glm.linear
done


module unload R/4.2.3


