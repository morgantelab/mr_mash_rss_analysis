#!/bin/bash

#SBATCH --job-name=run_merge_plink_mtag_output_ot_ukb
#SBATCH --cpus-per-task=2
#SBATCH --partition=fm-bigmem-1,fm-bigmem-2,fm-bigmem-3,fm-bigmem-4
#SBATCH --time=02:00:00
#SBATCH --mem=15G
#SBATCH --output=run_merge_plink_mtag_output_ot_ukb.%j.out
#SBATCH --error=run_merge_plink_mtag_output_ot_ukb.%j.err

###Load module
source /opt/intel/oneapi/mkl/2023.0.0/env/vars.sh intel64
module load R/4.2.3

###Set variables
export MKL_NUM_THREADS=1

TRAIT_NUM=0

for TRAIT in weight waist hip BMI TFM BFP DPa SPa;
do
    TRAIT_NUM=$[$TRAIT_NUM +1]

    Rscript ../code/merge_plink_mtag_output_ot_ukb.R --mtag_sumstats ../output/GWAS_for_regions/ukb_ot_chrAll_biallelic_snps_only_recoded_ids_info_06_gwas_for_regions.rsid.mtag_trait_"${TRAIT_NUM}".txt \
                                                     --plink_sumstats ../output/GWAS_for_regions/ukb_ot_chrAll_biallelic_snps_only_recoded_ids_info_06_gwas_for_regions.rsid."${TRAIT}".glm.linear \
           		                                       --ncores 2 \
               		                                   --output ../output/GWAS_for_regions/ukb_ot_chrAll_biallelic_snps_only_recoded_ids_info_06_gwas_for_regions_rsid_"${TRAIT}"_plink_mtag_combined.txt
done


module unload R/4.2.3


