#!/bin/bash

ml BCFtools/1.14-GCC-11.2.0

### Run script in cseqtl/results as sbatch slurm_geno.sh

# 1) Concatenate chromosomes into bcf file
# 2) Split into sample (indv) specific bcf files
# 3) Convert to tab delim format for aseq

GENO_DIR=/fh/scratch/delete90/kooperberg_c/topmed_freeze10/phased/whi_only
OUT_DIR=/fh/scratch/delete90/kooperberg_c/mjohnson/cseqtl/results/genotype

INFILE=freeze10b.whi_only_chr
# OUTFILE=WHI

# OUT_DIR2=${OUT_DIR}/split
# OUT_DIR3=${OUT_DIR2}/sub

# 1) Concat --------------------------------------------------------------------

# Filter particpants across chromsomes
cd $OUT_DIR

for sample in ${GENO_DIR}/${INFILE}{1..22}.bcf do;
do bcftools view --samples-file ${OUT_DIR}/geno_files_filt.txt -Oz ${OUT_DIR}/${sample};
done
