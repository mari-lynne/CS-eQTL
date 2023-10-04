#!/bin/bash

# Aims:
# Filter participants from LLS study from the larger WHI cohort WGS data

ml BCFtools/1.14-GCC-11.2.0

# 1) Set Directories and variables ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

GENO_DIR=/fh/scratch/delete90/kooperberg_c/topmed_freeze10/phased/whi_only # BCF File Input Directory
GENO_FILE=freeze10b.whi_only_chr

IN_DIR=/fh/scratch/delete90/kooperberg_c/mjohnson/cseqtl/results/genotype # Sample_file/manifest Directory
SAMPLE_FILE=lls_geno_filt.txt

OUT_DIR=${IN_DIR}/LLS/Sept

# Check if directories exist
if [ ! -d "$GENO_DIR" ] || [ ! -d "$IN_DIR" ] || [ ! -d "$OUT_DIR" ]; then
  echo "One or more directories do not exist. Exiting."
  exit 1
fi

# 2) Filter particpants across chromsomes ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Make a 1-dim array to store chromsome values (@ notation subsets array)
# For each chromsome submit a new job (23 jobs total)
# Then for each chromsome bcf file, subset participants and make a new bcf file
# -m2 -M2 -v snps selects biallelic snps

chr=( 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X )
for f in "${chr[@]}"
do
  sbatch -c 1 -t 3-0 --mem 20000 \
    --mail-user=mjohnso5@fredhutch.org \
    --mail-type=BEGIN,END,FAIL \
    --output=Rout/%j.out \
    --error=Rerr/%j.err \
    --wrap="bcftools view -Ob -m2 -M2 -v snps --force-samples \
            -S ${IN_DIR}/${SAMPLE_FILE} \
            ${GENO_DIR}/${GENO_FILE}${f}.bcf > \
            ${OUT_DIR}/${GENO_FILE}${f}.bcf"
done

echo "Sample Filtering - DONE!"
