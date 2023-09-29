#!/bin/bash

#SBATCH --job-name=compute_LD_bc_ukb
#SBATCH --cpus-per-task=10
#SBATCH --partition=fm-bigmem-1,fm-bigmem-2,fm-bigmem-3
#SBATCH --time=120:00:00
#SBATCH --mem=100G
#SBATCH --output=compute_LD_bc_ukb.%j.out
#SBATCH --error=compute_LD_bc_ukb.%j.err

###Get fold from arguments
FOLD=${1}

###Run ldstore
/data2/morgante_lab/fabiom/software/ldstore_v2.0_x86_64/ldstore_v2.0_x86_64 \
                                   --in-files ldstore_master_file_"${FOLD}" \
                                   --n-threads 10 \
                                   --write-text \
                                   --write-bdose \
                                   --bdose-version 1.0

