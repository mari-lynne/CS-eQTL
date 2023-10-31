#!/bin/bash
#SBATCH --job-name=cibersort_lls
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=20G
#SBATCH --output=Rout/cibersort_lls%_J.out
#SBATCH --output=Rerr/cibersort_lls%_J.err
#SBATCH --mail-type=begin,end,fail
#SBATCH --mail-user=mjohnso5@fredhutch.org

ml fhR/4.2.2-foss-2021b

Rscript 3_cibersort.R
