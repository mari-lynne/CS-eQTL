#!/bin/bash

### run in cseqtl/results as sbatch --test-only slurm_geno.sh

# 1) Concatenate chromosomes into bcf file
# 2) Split into sample (indv) specific bcf files
# 3) Convert to tab delim

GENO_DIR=/fh/scratch/delete90/kooperberg_c/topmed_freeze10/phased/whi_only
OUT_DIR=/fh/scratch/delete90/kooperberg_c/mjohnson/cseqtl/results/genotype
OUT_DIR2=${OUT_DIR}/split
FILENAME=WHI

# 1) Concat --------------------------------------------------------------------

bcftools concat ${OUT_DIR}/freeze10b.whi_only_chr{1..22}.bcf --naive -Ob -o ${OUT_DIR}/${FILENAME}_concat.bcf

# 2) Split ---------------------------------------------------------------------

bcftools +split ${OUT_DIR}/${FILENAME}_concat.bcf -Ob -o ${OUT_DIR2}

# convert to tab-delim
