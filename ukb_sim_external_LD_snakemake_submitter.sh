#!/bin/bash
#
#SBATCH --job-name=ukb_sim_sm
#SBATCH --ntasks=1
#SBATCH --partition=compute,bigmem
#SBATCH --time=336:00:00
#SBATCH --mem=2gb
#SBATCH --output=run/log/ukb_caucasian_white_british_unrel_100000_external_LD_equal_effects_indep_resid.%j.out
#SBATCH --error=run/log/ukb_caucasian_white_british_unrel_100000_external_LD_equal_effects_indep_resid.%j.err
#SBATCH --mail-type=fail
#SBATCH --mail-user=fabiom@clemson.edu

source /data2/morgante_lab/fabiom/software/miniconda3/etc/profile.d/conda.sh
source /data2/morgante_lab/fabiom/software/miniconda3/etc/profile.d/mamba.sh
conda activate mr_mash_rss_proj

#--dag | display | dot
#-p -n \
## test dag generation
# snakemake -p -n -s ukb_sim_external_LD_snakefile \
#           --configfile ukb_sim_external_LD_equal_effects_indep_resid.yaml
           #--rerun-triggers mtime \
           #--configfile ukb_sim_external_LD_equal_effects_50000causal_indep_resid.yaml 
           #--configfile ukb_sim_external_LD_trait_1_only_effects_indep_resid.yaml
           #--configfile ukb_sim_external_LD_blocks_shared_effects_indep_resid.yaml
           #--configfile ukb_sim_external_LD_equal_effects_10traits_indep_resid.yaml
           #--configfile ukb_sim_external_LD_equal_effects_low_pve_indep_resid.yaml

snakemake \
  -s ukb_sim_external_LD_snakefile \
  --profile slurm \
  --latency-wait 120 \
  -k \
  --configfile ukb_sim_external_LD_equal_effects_indep_resid.yaml
  #--configfile ukb_sim_external_LD_equal_effects_low_pve_indep_resid.yaml
  #--configfile ukb_sim_external_LD_blocks_shared_effects_indep_resid.yaml
  #--configfile ukb_sim_external_LD_trait_1_only_effects_indep_resid.yaml
  #--configfile ukb_sim_external_LD_equal_effects_50000causal_indep_resid.yaml
  #--configfile ukb_sim_external_LD_equal_effects_10traits_indep_resid.yaml
