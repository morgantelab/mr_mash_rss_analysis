#!/bin/bash

#SBATCH --job-name=run_obtain_finemapping_regions_ot_ukb
#SBATCH --cpus-per-task=4
#SBATCH --partition=fm-bigmem-1,fm-bigmem-2,fm-bigmem-3,fm-bigmem-4
#SBATCH --time=02:00:00
#SBATCH --mem=10G
#SBATCH --output=run_obtain_finemapping_regions_ot_ukb.%j.out
#SBATCH --error=run_obtain_finemapping_regions_ot_ukb.%j.err

###Get chromosome number from array task id
CHR="${SLURM_ARRAY_TASK_ID}"

###Load module
source /opt/intel/oneapi/mkl/2023.0.0/env/vars.sh intel64
module load R/4.2.3

###Set variables
export MKL_NUM_THREADS=1

Rscript ../code/obtain_finemapping_regions_ot_ukb.R --sumstats_prefix ../output/GWAS_for_regions/ukb_ot_chrAll_biallelic_snps_only_recoded_ids_info_06_gwas_for_regions_rsid_ \
                                                    --sumstats_suffix _plink_mtag_combined.txt \
                                                    --chr "${CHR}" \
                                                    --traits "c('BMI','BMR','DPa','SPa','WHR','height','hip')" \
                                                    --region_size 250000 \
                                                    --mhc FALSE \
                                                    --sig_threshold 5e-8 \
                                                    --gwas_method mtag \
                                                    --ncores 4 \
                                                    --output_prefix_snp ../data/other_traits_finemapping_regions/othertraits \
                                                    --output_prefix_boundaries ../data/misc/othertraits


module unload R/4.2.3


