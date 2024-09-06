#!/bin/bash

#SBATCH --job-name=run_mtag_ot_ukb
#SBATCH --cpus-per-task=8
#SBATCH --partition=fm-bigmem-1,fm-bigmem-2,fm-bigmem-3,fm-bigmem-4
#SBATCH --time=02:00:00
#SBATCH --mem=100G
#SBATCH --output=run_mtag_ot_ukb.%j.out
#SBATCH --error=run_mtag_ot_ukb..%j.err

###Load module
module load python2.7-fm/2.7.18

###Run MTAG
python2 /data2/morgante_lab/fabiom/software/mtag/mtag.py  \
    --sumstats ../output/GWAS_for_regions/ukb_ot_chrAll_biallelic_snps_only_recoded_ids_info_06_gwas_for_regions.rsid.weight.glm.linear,../output/GWAS_for_regions/ukb_ot_chrAll_biallelic_snps_only_recoded_ids_info_06_gwas_for_regions.rsid.waist.glm.linear,../output/GWAS_for_regions/ukb_ot_chrAll_biallelic_snps_only_recoded_ids_info_06_gwas_for_regions.rsid.hip.glm.linear,../output/GWAS_for_regions/ukb_ot_chrAll_biallelic_snps_only_recoded_ids_info_06_gwas_for_regions.rsid.BMI.glm.linear,../output/GWAS_for_regions/ukb_ot_chrAll_biallelic_snps_only_recoded_ids_info_06_gwas_for_regions.rsid.TFM.glm.linear,../output/GWAS_for_regions/ukb_ot_chrAll_biallelic_snps_only_recoded_ids_info_06_gwas_for_regions.rsid.BFP.glm.linear,../output/GWAS_for_regions/ukb_ot_chrAll_biallelic_snps_only_recoded_ids_info_06_gwas_for_regions.rsid.DPa.glm.linear,../output/GWAS_for_regions/ukb_ot_chrAll_biallelic_snps_only_recoded_ids_info_06_gwas_for_regions.rsid.SPa.glm.linear \
    --snp_name rsID \
    --chr_name \#CHROM \
    --bpos_name POS \
    --a1_name A1 \
    --a2_name OMITTED \
    --eaf_name A1_FREQ \
    --z_name T_STAT \
    --beta_name BETA \
    --se_name SE \
    --n_name OBS_CT \
    --p_name P \
    --maf_min 0 \
    --n_min 0 \
    --out ../output/GWAS_for_regions/ukb_ot_chrAll_biallelic_snps_only_recoded_ids_info_06_gwas_for_regions.rsid.mtag

###Unload module
module unload python2.7-fm/2.7.18



