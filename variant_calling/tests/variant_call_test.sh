#!/bin/bash

ml SAMtools
ml BCFtools
ml GATK
ml picard

# Define the reference genome file
REFERENCE="/fh/scratch/delete90/kooperberg_c/mjohnson/cseqtl/topmed_variant_calling/test_data/resources/ref/hs38DH.fa"

# Define the output directory
OUT_DIR="/fh/scratch/delete90/kooperberg_c/mjohnson/cseqtl/topmed_variant_calling/test_data"

# Define the input directory with the CRAM files
IN_DIR="/fh/scratch/delete90/kooperberg_c/mjohnson/cseqtl/topmed_variant_calling/test_data"

# Create an array with the names of the CRAM files
FILES=($(ls $IN_DIR/*.cram))

for FILE in "${FILES[@]}"; do
# Extract the sample name from the filename
SAMPLE_NAME=$(basename "$FILE" .src.cram)

# Define the output VCF file
VCF="$OUT_DIR/$SAMPLE_NAME.vcf"

# Create a pileup and call variants with bcftools
samtools mpileup -uf $REFERENCE $FILE | bcftools call -mv -Ov -o $VCF

# Validate the VCF file with Picard
java -jar picard.jar ValidateVariants \
I=$VCF \
R=$REFERENCE \
IGNORE=ErrorType

done