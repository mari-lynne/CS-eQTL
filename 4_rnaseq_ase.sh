#!/bin/bash

#SBATCH --array=[1-21]
#SBATCH --nodes=1
#SBATCH --mem-per-cpu=32000
#SBATCH --cpus-per-task=1
#SBATCH --output=Rout/par-%j-%a.out
#SBATCH --error=Rout/par-%j-%a.err
#SBATCH --mail-type=begin
#SBATCH --mail-type=end
#SBATCH --mail-type=fail
#SBATCH --mail-user=mjohnso5@fredhutch.org

ml fhR/4.2.2.1-foss-2021b

Rscript 4_rnaseq_ase.R $SLURM_JOB_NAME

# Run in fh/scratch/delete90/kooperberg_c/mjohnson/cseqtl as sbatch -J manifest_rna_geno.csv 4_rnaseq_ase.sh

