#!/bin/bash

#SBATCH --job-name=compute_number_of_variants_per_chromosome
#SBATCH --cpus-per-task=1
#SBATCH --partition=fm-bigmem-1,fm-bigmem-2,fm-bigmem-3
#SBATCH --time=02:00:00
#SBATCH --mem=10G
#SBATCH --output=compute_number_of_variants_per_chromosome.%j.out
#SBATCH --error=compute_number_of_variants_per_chromosome.%j.err


###Compute allele frequencies per chromosome
for i in {1..22}
do
   /data2/morgante_lab/fabiom/software/plink19 --bfile ../data/genotypes/ukb_caucasian_white_british_unrel_100000 \
					       						--keep-allele-order \
                           						--chr ${i} \
                           						--freq \
                           						--out ../data/LD_matrices/number_of_variants/ukb_caucasian_white_british_unrel_100000_chr"${i}" \
                           						--threads 1 
done

###Load module
source /opt/intel/oneapi/mkl/2023.0.0/env/vars.sh intel64
module load R/4.2.3

###Run script to compute number of variants per chromosome
Rscript ../code/compute_num_variants_per_chr_ukb.R --input ukb_caucasian_white_british_unrel_100000

###Unload module
module unload R/4.2.3

###Remove intermediate files
rm -rf ../data/LD_matrices/number_of_variants/*.frq ../data/LD_matrices/number_of_variants/*.log
