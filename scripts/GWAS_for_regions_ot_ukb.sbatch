#!/bin/bash

#SBATCH --job-name=GWAS_for_regions_ot_ukb
#SBATCH --cpus-per-task=8
#SBATCH --partition=fm-bigmem-1,fm-bigmem-2,fm-bigmem-3,fm-bigmem-4
#SBATCH --time=12:00:00
#SBATCH --mem=310G
#SBATCH --output=GWAS_for_regions_ot_ukb.%j.out
#SBATCH --error=GWAS_for_regions_ot_ukb.%j.err

###Get chromosome number from array task id
CHR="${SLURM_ARRAY_TASK_ID}"

###Load module
source /opt/intel/oneapi/mkl/2023.0.0/env/vars.sh intel64
module load R/4.2.3

###Set variables
export MKL_NUM_THREADS=1
export OMP_NUM_THREADS=1


Rscript ../code/filter_variants_by_info_score_ukb.R --mfi /data2/morgante_lab/data/ukbiobank/genotypes/imputed/ukb_mfi_chr"${CHR}"_v3.txt \
                                                    --info_thresh 0.6 \
                                                    --chr ${CHR} \
                                                    --ncores 8 \
                                                    --output ../data/misc/ukb_imp_chr"${CHR}"_variants_info_06.txt

###Keep only biallelic SNPs
/data2/morgante_lab/fabiom/software/plink2 \
    --bgen /data2/morgante_lab/data/ukbiobank/genotypes/imputed/ukb22828_c"${CHR}"_b0_v3.bgen ref-first \
    --sample /data2/morgante_lab/data/ukbiobank/genotypes/imputed/ukb22828_c"${CHR}"_b0_v3_s487271.sample \
    --snps-only \
    --max-alleles 2 \
    --make-pgen \
    --out ../data/genotypes/ukb22828_c"${CHR}"_b0_v3_biallelic_snps_only \
    --memory 310000 \
    --threads 8
    
###Recode SNP ids so that they are unique and write a new PVAR file    
/data2/morgante_lab/fabiom/software/plink2 \
    --pfile ../data/genotypes/ukb22828_c"${CHR}"_b0_v3_biallelic_snps_only \
    --set-all-var-ids @:#_\$r_\$a \
    --make-just-pvar \
    --out ../data/genotypes/ukb22828_c"${CHR}"_b0_v3_biallelic_snps_only_recoded_ids \
    --memory 310000 \
    --threads 8

###Perform GWAS for variants with INFO >= 0.6, MAF > 0.001
/data2/morgante_lab/fabiom/software/plink2 \
    --pgen ../data/genotypes/ukb22828_c"${CHR}"_b0_v3_biallelic_snps_only.pgen \
    --pvar ../data/genotypes/ukb22828_c"${CHR}"_b0_v3_biallelic_snps_only_recoded_ids.pvar \
    --psam ../data/genotypes/ukb22828_c"${CHR}"_b0_v3_biallelic_snps_only.psam \
    --extract ../data/misc/ukb_imp_chr"${CHR}"_variants_info_06.txt \
    --keep ../data/misc/ukb_cleaned_ot_ind_ids.txt \
    --pheno ../data/phenotypes/ukb_cleaned_ot_adjusted_pheno_plink.txt \
    --glm no-x-sex omit-ref allow-no-covars \
    --maf 0.001 \
    --out ../output/GWAS_for_regions/ukb_ot_chr"${CHR}"_biallelic_snps_only_recoded_ids_info_06_gwas_for_regions \
    --memory 310000 \
    --threads 8

module unload R/4.2.3


