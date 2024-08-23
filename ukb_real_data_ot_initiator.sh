#!/bin/bash
mkdir -p ./run/{log,logs_slurm} | sbatch ./ukb_real_data_ot_snakemake_submitter.sh
