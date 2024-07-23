#!/bin/bash

#SBATCH --job-name=process_1KG_data_for_LD
#SBATCH --cpus-per-task=1
#SBATCH --partition=fm-bigmem-1,fm-bigmem-2,fm-bigmem-3,fm-bigmem-4,compute,compute,bigmem
#SBATCH --time=24:00:00
#SBATCH --mem=10G
#SBATCH --output=process_1KG_data_for_LD.%j.out
#SBATCH --error=process_1KG_data_for_LD.%j.err

###Get chromosome number from array task id
CHR="${SLURM_ARRAY_TASK_ID}"

###Load modules
source /opt/intel/oneapi/mkl/2023.0.0/env/vars.sh intel64
module load R/4.2.3
module load bcftools/1.9

###Get variants info from UKB LD data
cut -f1,4-6 "../data/genotypes/ukb_caucasian_white_british_unrel_for_LD_genetic_dist_chr${CHR}.bim" > "/scratch3/fabiom/tmp/ukb_caucasian_white_british_unrel_for_LD_genetic_dist_variants_chr${CHR}.txt"

###Filter VCF for variants in common with ukb_caucasian_white_british_unrel_for_LD_genetic_dist based on POS REF ALT (https://www.biostars.org/p/170965/)
bcftools view -Ov -R "/scratch3/fabiom/tmp/ukb_caucasian_white_british_unrel_for_LD_genetic_dist_variants_chr${CHR}.txt" "/data2/morgante_lab/data/1000G/ALL.chr${CHR}.phase3_integrated.20130502.genotypes.bcf" \
 | grep -Ef <(awk 'BEGIN{FS=OFS="\t";print "#"};{print "^"$1,$2,"[^\t]+",$3,$4"\t"}' "/scratch3/fabiom/tmp/ukb_caucasian_white_british_unrel_for_LD_genetic_dist_variants_chr${CHR}.txt") \
 > "/scratch3/fabiom/tmp/ALL.chr${CHR}.phase3_integrated.20130502.genotypes.filtered.vcf"
 
###Annotate with variant ID from UKB
Rscript ../code/prepare_variants_info_file_ukb.R --input "../data/genotypes/ukb_caucasian_white_british_unrel_for_LD_genetic_dist_chr${CHR}.bim" \
                                                 --output "../data/genotypes/ukb_caucasian_white_british_unrel_for_LD_genetic_dist_chr${CHR}_variants_info.txt"
                                         
bgzip "../data/genotypes/ukb_caucasian_white_british_unrel_for_LD_genetic_dist_chr${CHR}_variants_info.txt"
tabix -s1 -b2 -e2 "../data/genotypes/ukb_caucasian_white_british_unrel_for_LD_genetic_dist_chr${CHR}_variants_info.txt.gz"

bcftools annotate -a "../data/genotypes/ukb_caucasian_white_british_unrel_for_LD_genetic_dist_chr${CHR}_variants_info.txt.gz" -c CHROM,POS,REF,ALT,ID \
  -Ov -o "/scratch3/fabiom/tmp/ALL.chr${CHR}.phase3_integrated.20130502.genotypes.filtered.ids.vcf" "/scratch3/fabiom/tmp/ALL.chr${CHR}.phase3_integrated.20130502.genotypes.filtered.vcf"

###Transform VCF to plink format and only keep EUR individuals
/data2/morgante_lab/fabiom/software/plink2 \
	--vcf "/scratch3/fabiom/tmp/ALL.chr${CHR}.phase3_integrated.20130502.genotypes.filtered.ids.vcf" \
  --ref-allele force "/scratch3/fabiom/tmp/ALL.chr${CHR}.phase3_integrated.20130502.genotypes.filtered.ids.vcf" 5 3 '#' \
	--keep /data2/morgante_lab/data/1000G/bigsnpr/EUR_ids.txt \
	--make-bed \
	--out "/scratch3/fabiom/tmp/1000G_phase3_eur_for_LD_chr${CHR}"

###AddGenetic distance to plink file	
/data2/morgante_lab/fabiom/software/plink19 --bfile "/scratch3/fabiom/tmp/1000G_phase3_eur_for_LD_chr${CHR}" \
                                            --keep-allele-order \
                                            --cm-map ../data/misc/1000-genomes-genetic-map-ldak/genetic_map_b37/genetic_map_chr@_combined_b37.txt \
                                            --make-bed \
                                            --out "../data/genotypes/1000G_phase3_eur_for_LD_genetic_dist_chr${CHR}"

module unload R/4.2.3
module unload bcftools/1.9

###Delete intermediate files
rm "/scratch3/fabiom/tmp/1000G_phase3_eur_for_LD_chr${CHR}.bed" "/scratch3/fabiom/tmp/1000G_phase3_eur_for_LD_chr${CHR}.bim" "/scratch3/fabiom/tmp/1000G_phase3_eur_for_LD_chr${CHR}.fam" \
"/scratch3/fabiom/tmp/ukb_caucasian_white_british_unrel_for_LD_genetic_dist_variants_chr${CHR}.txt" "/scratch3/fabiom/tmp/1000G_phase3_eur_for_LD_chr${CHR}.log" \
"../data/genotypes/1000G_phase3_eur_for_LD_genetic_dist_chr${CHR}.log" 

