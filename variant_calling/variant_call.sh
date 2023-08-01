#!/bin/bash

# Define the reference genome file
REFERENCE="/path/to/reference/genome.fasta"

# Define the output directory
OUT_DIR="/path/to/output/directory/"

# Define the input directory with the CRAM files
IN_DIR="/fh/scratch/delete90/kooperberg_c/wgs_bam"

# Create an array with the names of the CRAM files
FILES=($(ls $IN_DIR/*.cram))

for FILE in "${FILES[@]}"; do
    # Extract the sample name from the filename
    SAMPLE_NAME=$(basename "$FILE" .src.cram)

    # Define the output VCF file
    VCF="$OUT_DIR/$SAMPLE_NAME.vcf"

    # Create a pileup and call variants with bcftools
    samtools mpileup -uf $REFERENCE $FILE | bcftools call -mv -Ov -o $VCF

    # Alternatively, use GATK HaplotypeCaller for variant calling (uncomment the below line if you want to use GATK instead of bcftools)
    # java -jar $GATK -T HaplotypeCaller -R $REFERENCE -I $FILE -o $VCF

    # Validate the VCF file with Picard
    java -jar picard.jar ValidateVariants \
         I=$VCF \
         R=$REFERENCE \
         IGNORE=ErrorType

done