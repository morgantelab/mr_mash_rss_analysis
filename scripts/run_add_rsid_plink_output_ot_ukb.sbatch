#!/bin/bash

#SBATCH --job-name=run_add_rsid_plink_output_ot_ukb
#SBATCH --cpus-per-task=8
#SBATCH --partition=fm-bigmem-1,fm-bigmem-2,fm-bigmem-3,fm-bigmem-4
#SBATCH --time=02:00:00
#SBATCH --mem=10G
#SBATCH --output=run_add_rsid_plink_output_ot_ukb.%j.out
#SBATCH --error=run_add_rsid_plink_output_ot_ukb..%j.err

###Get chromosome number from array task id
CHR="${SLURM_ARRAY_TASK_ID}"

###Load module
source /opt/intel/oneapi/mkl/2023.0.0/env/vars.sh intel64
module load R/4.2.3

###Set variables
export MKL_NUM_THREADS=1

for TRAIT in BMI BMR WHR height hip DPa SPa;
do
    Rscript ../code/add_rsid_plink_output_ot_ukb.R --mfi /data2/morgante_lab/data/ukbiobank/genotypes/imputed/ukb_mfi_chr"${CHR}"_v3.txt \
                                                   --plink_sumstats ../output/GWAS_for_regions/ukb_ot_chr"${CHR}"_biallelic_snps_only_recoded_ids_info_06_gwas_for_regions."${TRAIT}".glm.linear \
       		                                   --chr "${CHR}" \
           		                           --ncores 8 \
               		                           --output ../output/GWAS_for_regions/ukb_ot_chr"${CHR}"_biallelic_snps_only_recoded_ids_info_06_gwas_for_regions.rsid."${TRAIT}".glm.linear
done


module unload R/4.2.3


