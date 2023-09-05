#!/bin/bash

# Steps ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# 1) Filter WHI particpants across chromsomes

ml BCFtools/1.14-GCC-11.2.0

# Set Directories and variables ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Input
GENO_DIR=/fh/scratch/delete90/kooperberg_c/topmed_freeze10/phased/whi_only
# Where the sample_file is saved
IN_DIR=/fh/scratch/delete90/kooperberg_c/mjohnson/cseqtl/results/genotype
SAMPLE_FILE=lls_geno_filt.txt
GENO_FILE=freeze10b.whi_only_chr

# OUTPUT
OUT_DIR=${IN_DIR}/LLS

# Start Script
# Check if directories exist
if [ ! -d "$GENO_DIR" ] || [ ! -d "$IN_DIR" ] || [ ! -d "$OUT_DIR" ]; then
  echo "One or more directories do not exist. Exiting."
  exit 1
fi

# 1) Filter particpants across chromsomes ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# For each chromsome bcf file subset participants an make a new bcf file
# Each chromsome submits an new job (so we will have 23 jobs total)
# Use @ notation to subset 1 dim array

chr=( 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X )
for f in "${chr[@]}"
do
  sbatch -c 1 -t 3-0 --mem 20000 \
    --mail-user=mjohnso5@fredhutch.org \
    --mail-type=BEGIN,END,FAIL \
    --output=Rout/%j.out \
    --error=Rerr/%j.err \
    --wrap="bcftools view -Ob --force-samples \
            -S ${IN_DIR}/${SAMPLE_FILE} \
            ${GENO_DIR}/${GENO_FILE}${f}.bcf > \
            ${OUT_DIR}/${GENO_FILE}${f}.bcf"
done

echo "Sample Filtering DONE"
