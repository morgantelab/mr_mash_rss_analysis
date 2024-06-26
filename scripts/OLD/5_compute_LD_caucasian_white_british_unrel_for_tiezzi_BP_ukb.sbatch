#!/bin/bash

#SBATCH --job-name=compute_LD_caucasian_white_british_unrel_for_tiezzi_BP
#SBATCH --cpus-per-task=16
#SBATCH --partition=fm-bigmem-1,fm-bigmem-2,fm-bigmem-3
#SBATCH --time=05:00:00
#SBATCH --mem=210G
#SBATCH --output=compute_LD_caucasian_white_british_unrel_for_tiezzi_BP.%j.out
#SBATCH --error=compute_LD_caucasian_white_british_unrel_for_tiezzi_BP.%j.err

###Get chromosome number from array task id
CHR="${SLURM_ARRAY_TASK_ID}"

###Run plink
echo "Start analyzing chromosome ${CHR} ..."

###Create temporary filtered file in BGEN v1.1
/data2/morgante_lab/fabiom/software/plink2 --bgen /data2/morgante_lab/data/ukbiobank/genotypes/imputed/ukb22828_c"${CHR}"_b0_v3.bgen ref-first \
                                           --sample /data2/morgante_lab/data/ukbiobank/genotypes/imputed/ukb22828_c"${CHR}"_b0_v3_s487271.sample \
                                           --keep ../data/LD_matrices/ukb_tiezzi_BP_inds_for_LD.txt \
                                           --extract ../data/genotypes/ukb_hm3_variant_rsids.txt \
                                           --export bgen-1.1 \
                                           --out /scratch1/fabiom/tmp/temp_data_"${CHR}" \
                                           --threads 16

###Compute LD matrix
/data2/morgante_lab/fabiom/software/plink19 --bgen /scratch1/fabiom/tmp/temp_data_"${CHR}".bgen \
                                            --sample /scratch1/fabiom/tmp/temp_data_"${CHR}".sample \
                                            --keep-allele-order \
                                            --r bin 'yes-really' \
                                            --out ../data/LD_matrices/ukb_caucasian_white_british_unrel_for_tiezzi_BP_chr"${CHR}"_LD \
                                            --threads 16

###Remove intermediate files
rm /scratch1/fabiom/tmp/temp_data_"${CHR}".*

echo "Finished analyzing chromosome ${CHR}."
