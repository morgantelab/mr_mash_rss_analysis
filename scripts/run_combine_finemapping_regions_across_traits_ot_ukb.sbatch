#!/bin/bash

#SBATCH --job-name=run_combine_finemapping_regions_across_traits_ot_ukb
#SBATCH --cpus-per-task=4
#SBATCH --partition=fm-bigmem-1,fm-bigmem-2,fm-bigmem-3,fm-bigmem-4
#SBATCH --time=02:00:00
#SBATCH --mem=15G
#SBATCH --output=run_combine_finemapping_regions_across_traits_ot_ukb.%j.out
#SBATCH --error=run_combine_finemapping_regions_across_traits_ot_ukb.%j.err

###Load module
source /opt/intel/oneapi/mkl/2023.0.0/env/vars.sh intel64
module load R/4.2.3

###Set variables
export MKL_NUM_THREADS=1

Rscript ../code/combine_finemapping_regions_across_traits_ot_ukb.R --regions_prefix ../data/misc/othertraits_ \
                                                                  --regions_suffix _regions.txt \
                                                                  --traits "c('weight', 'waist', 'hip', 'BMI', 'TFM', 'BFP', 'DPa', 'SPa')" \
                                                                  --ncores 4 \
                                                                  --output_prefix_snp ../data/other_traits_finemapping_regions/othertraits \
                                                                  --output_boundaries ../data/misc/othertraits_combined_regions_chrAll_info.txt


module unload R/4.2.3


