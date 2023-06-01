#!/bin/bash
mkdir -p /data2/morgante_lab/ukbiobank_projects/mr_mash_rss/run/{log,logs_slurm} | sbatch ./ukb_sim_snakemake_submitter.sh
