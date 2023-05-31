#!/bin/bash
#
#SBATCH --job-name=sm_test
#SBATCH --ntasks=1
#SBATCH --partition=compute
#SBATCH --time=72:00:00
#SBATCH --mem=2gb
#SBATCH --output=run/log/mr_mash_rss_pipeline_sample_filter_ukb_caucasian_white_british_unrel.%j.out
#SBATCH --error=run/log/mr_mash_rss_pipeline_sample_filter_ukb_caucasian_white_british_unrel.%j.err
#SBATCH --mail-type=fail
#SBATCH --mail-user=fabiom@clemson.edu

source /data2/morgante_lab/fabiom/software/miniconda3/etc/profile.d/conda.sh
source /data2/morgante_lab/fabiom/software/miniconda3/etc/profile.d/mamba.sh
conda activate mr_mash_rss_proj

#--dag | display | dot
#-p -n \
## test dag generation
snakemake -p -n -s Snakefile --configfile ukb_mrmash.yaml

# snakemake \
# -s Snakefile \
# --profile slurm \
# --configfile ukb_mrmash.yaml \
# --latency-wait 120
