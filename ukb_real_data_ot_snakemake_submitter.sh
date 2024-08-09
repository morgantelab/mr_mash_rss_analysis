#!/bin/bash
#
#SBATCH --job-name=ukb_real_data_ot_sm
#SBATCH --ntasks=1
#SBATCH --partition=compute
#SBATCH --time=336:00:00
#SBATCH --mem=2gb
#SBATCH --output=run/log/ukb_real_data_ot.%j.out
#SBATCH --error=run/log/ukb_real_data_ot.%j.err

source /data2/morgante_lab/fabiom/software/miniconda3/etc/profile.d/conda.sh
source /data2/morgante_lab/fabiom/software/miniconda3/etc/profile.d/mamba.sh
conda activate mr_mash_rss_proj

#--dag | display | dot
#-p -n \
## test dag generation
#snakemake -p -n -s ukb_real_data_ot_snakefile \
#         --configfile ukb_real_data_ot.yaml
           #--rerun-triggers mtime \

snakemake \
  -s ukb_real_data_ot_snakefile \
  --profile slurm \
  --latency-wait 120 \
  --configfile ukb_real_data_ot.yaml
