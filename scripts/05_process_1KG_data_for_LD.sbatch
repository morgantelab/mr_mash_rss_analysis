#!/bin/bash

#SBATCH --job-name=process_1KG_data_for_LD
#SBATCH --cpus-per-task=1
#SBATCH --partition=fm-bigmem-1,fm-bigmem-2,fm-bigmem-3,fm-bigmem-4,compute,compute,bigmem
#SBATCH --time=48:00:00
#SBATCH --mem=10G
#SBATCH --output=process_1KG_data_for_LD.%j.out
#SBATCH --error=process_1KG_data_for_LD.%j.err

###Get chromosome number from array task id
CHR="${SLURM_ARRAY_TASK_ID}"

###Load modules
source /opt/intel/oneapi/mkl/2023.0.0/env/vars.sh intel64
module load R/4.2.3
module load bcftools/1.9

###Create variants info VCF
Rscript ../code/prepare_variants_info_file_ukb.R --input "../data/genotypes/ukb_caucasian_white_british_unrel_for_LD_genetic_dist_chr${CHR}.bim" \
                                                 --output "/scratch3/fabiom/tmp/ukb_caucasian_white_british_unrel_for_LD_genetic_dist_chr${CHR}_variants_info.txt"
                                                 
echo "##fileformat=VCFv4.3" | cat - "/scratch3/fabiom/tmp/ukb_caucasian_white_british_unrel_for_LD_genetic_dist_chr${CHR}_variants_info.txt" \
      > "/scratch3/fabiom/tmp/ukb_caucasian_white_british_unrel_for_LD_genetic_dist_chr${CHR}_variants_info.vcf"
                                         
bgzip "/scratch3/fabiom/tmp/ukb_caucasian_white_british_unrel_for_LD_genetic_dist_chr${CHR}_variants_info.vcf"
tabix -s1 -b2 -e2 "/scratch3/fabiom/tmp/ukb_caucasian_white_british_unrel_for_LD_genetic_dist_chr${CHR}_variants_info.vcf.gz"

###Filter BCF for variants in common with ukb_caucasian_white_british_unrel_for_LD_genetic_dist based on CHROM POS REF ALT 
bcftools isec -Oz -c none -p "/scratch3/fabiom/tmp/ALL.chr${CHR}.phase3_integrated.20130502.genotypes.filtered" \
              -n=2 -w2 "/scratch3/fabiom/tmp/ukb_caucasian_white_british_unrel_for_LD_genetic_dist_chr${CHR}_variants_info.vcf.gz" \
              "/data2/morgante_lab/data/1000G/ALL.chr${CHR}.phase3_integrated.20130502.genotypes.bcf"

###Annotate with variant ID from UKB
bcftools annotate -a "/scratch3/fabiom/tmp/ukb_caucasian_white_british_unrel_for_LD_genetic_dist_chr${CHR}_variants_info.vcf.gz" -c CHROM,POS,REF,ALT,ID \
  -Ov -o "/scratch3/fabiom/tmp/ALL.chr${CHR}.phase3_integrated.20130502.genotypes.filtered.ids.vcf" "/scratch3/fabiom/tmp/ALL.chr${CHR}.phase3_integrated.20130502.genotypes.filtered/0001.vcf.gz"

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
rm -rf "/scratch3/fabiom/tmp/1000G_phase3_eur_for_LD_chr${CHR}.bed" "/scratch3/fabiom/tmp/1000G_phase3_eur_for_LD_chr${CHR}.bim" "/scratch3/fabiom/tmp/1000G_phase3_eur_for_LD_chr${CHR}.fam" \
"/scratch3/fabiom/tmp/1000G_phase3_eur_for_LD_chr${CHR}.log" "ALL.chr${CHR}.phase3_integrated.20130502.genotypes.filtered.ids.vcf" \
"ukb_caucasian_white_british_unrel_for_LD_genetic_dist_chr${CHR}_variants_info.txt" "ukb_caucasian_white_british_unrel_for_LD_genetic_dist_chr${CHR}_variants_info.vcf.gz" \
"ukb_caucasian_white_british_unrel_for_LD_genetic_dist_chr${CHR}_variants_info.vcf.gz.tbi"
