#!/bin/bash

#SBATCH --job-name=run_obtain_finemapping_regions_single_trait_ot_ukb
#SBATCH --cpus-per-task=4
#SBATCH --partition=fm-bigmem-1,fm-bigmem-2,fm-bigmem-3,fm-bigmem-4
#SBATCH --time=02:00:00
#SBATCH --mem=15G
#SBATCH --output=run_obtain_finemapping_regions_single_trait_ot_ukb.%j.out
#SBATCH --error=run_obtain_finemapping_regions_single_trait_ot_ukb.%j.err

###Get chromosome number from array task id
TRAIT="${SLURM_ARRAY_TASK_ID}"

###Load module
source /opt/intel/oneapi/mkl/2023.0.0/env/vars.sh intel64
module load R/4.2.3

###Set variables
export MKL_NUM_THREADS=1

###
TRAITS=("weight" "waist" "hip" "BMI" "TFM" "BFP" "DPa" "SPa")

TRAIT_NAME="${TRAITS[${TRAIT}]}"

Rscript ../code/obtain_finemapping_regions_single_trait_ot_ukb.R --sumstats ../output/GWAS_for_regions/ukb_ot_chrAll_biallelic_snps_only_recoded_ids_info_06_gwas_for_regions_rsid_"${TRAIT_NAME}"_plink_mtag_combined.txt \
                                                                --region_size 250000 \
                                                                --mhc FALSE \
                                                                --sig_threshold 5e-8 \
                                                                --gwas_method mtag \
                                                                --ncores 4 \
                                                                --output ../data/misc/othertraits_"${TRAIT_NAME}"_regions.txt


module unload R/4.2.3


