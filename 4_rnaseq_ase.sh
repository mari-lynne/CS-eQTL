#!/bin/bash

#SBATCH --array=[1-1327]
#SBATCH --job-name=aseq_data_prep
#SBATCH --nodes=1
#SBATCH --mem-per-cpu=20G
#SBATCH --output=Rout/par-%j-%a.out
#SBATCH --error=Rerr/par-%j-%a.err
#SBATCH --mail-type=begin,end,fail
#SBATCH --mail-user=mjohnso5@fredhutch.org

ml fhR/4.2.2.1-foss-2021b

Rscript 2_rnaseq_ase.R $SLURM_JOB_NAME

# Subset = Run in fh/scratch/delete90/kooperberg_c/mjohnson/cseqtl as sbatch -J manifest_rna_geno.csv 4_rnaseq_ase.sh
# LLS Run in fh/scratch/delete90/kooperberg_c/mjohnson/cseqtl/scripts/rnaseq_scripts as sbatch -J lls_rna_array.csv 2_rnaseq_ase.sh

