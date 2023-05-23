#!/bin/bash

#SBATCH --nodes=1
#SBATCH --cpus-per-task=4
#SBATCH --output=Rout/par-%j-%a.out   # File to which STDOUT will be written, including job ID (%j)
#SBATCH --error=Rout/par-%j-%a.out    # File to which STDERR will be written, including job ID (%j)
#SBATCH --mail-type=begin
#SBATCH --mail-type=end
#SBATCH --mail-type=fail
#SBATCH --mail-user=mjohnso5@fredhutch.org

### run in cseqtl/results as sbatch --test-only -c 1 -t 2-0 --mem 8000 geno_prep.sh

# 1) Concatenate chromosomes into bcf file
# 2) Split into sample (indv) specific bcf files
# 3) Convert to tab delim

GENO_DIR=/fh/scratch/delete90/kooperberg_c/topmed_freeze10/phased/whi_only
OUT_DIR=/fh/scratch/delete90/kooperberg_c/mjohnson/cseqtl/results/genotype
OUT_DIR2=${OUT_DIR}/split
FILENAME=WHI

ml BCFtools/1.14-GCC-11.2.0

# 1) Concat --------------------------------------------------------------------

bcftools concat ${OUT_DIR}/freeze10b.whi_only_chr{1..22}.bcf --naive --write.index -Ob -o ${OUT_DIR}/${FILENAME}_concat.bcf

# 2) Split ---------------------------------------------------------------------

bcftools +split ${OUT_DIR}/${FILENAME}_concat.bcf -Ob -o ${OUT_DIR2}

# convert to tab-delim