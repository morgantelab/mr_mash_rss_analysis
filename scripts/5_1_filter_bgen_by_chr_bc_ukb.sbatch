#!/bin/bash

#SBATCH --job-name=filter_bgen_by_chr_bc_ukb
#SBATCH --cpus-per-task=1
#SBATCH --partition=fm-bigmem-1,fm-bigmem-2,fm-bigmem-3
#SBATCH --time=24:00:00
#SBATCH --mem=2G
#SBATCH --output=filter_bgen_by_chr_bc_ukb.%j.out
#SBATCH --error=filter_bgen_by_chr_bc_ukb.%j.err

###Get chromosome number from array task id
CHR="${SLURM_ARRAY_TASK_ID}"

###Get fold from arguments
FOLD=${1}

###Filter bgen
echo "Start analyzing chromosome ${CHR} and fold ${FOLD} ..."


/data2/morgante_lab/fabiom/software/qctool/build/release/apps/qctool_v2.2.0 \
                                   -g /data2/morgante_lab/data/ukbiobank/genotypes/imputed/ukb22828_c"${CHR}"_b0_v3.bgen \
                                   -s /data2/morgante_lab/data/ukbiobank/genotypes/imputed/ukb22828_c"${CHR}"_b0_v3_s487271.sample \
                                   -og ../tmp/ukb22828_c"${CHR}"_b0_v3_bc_filtered_"${FOLD}".bgen \
                                   -os ../tmp/ukb22828_c"${CHR}"_b0_v3_bc_filtered_"${FOLD}".sample \
                                   -incl-samples ../data/LDstore_input_files/ukb_bc_training_samples_LDstore_"${FOLD}".incl \
                                   -incl-rsids ../data/LDstore_input_files/ukb_bc_chr"${CHR}"_LDstore_rsid.txt

###Produce bgi
cd ../tmp
/data2/morgante_lab/fabiom/software/bgen/build/apps/bgenix -g ukb22828_c"${CHR}"_b0_v3_bc_filtered_"${FOLD}".bgen \
                                                          -index

cd ../run

echo "Finished analyzing chromosome ${CHR} and fold ${FOLD}."

