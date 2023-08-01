#!/bin/bash

# Define the reference genome file
REFERENCE="/path/to/reference/genome.fasta"

# Define the input directory with the VCF files
VCF_DIR="/path/to/output/directory/"

# Create an array with the names of the VCF files
FILES=($(ls $VCF_DIR/*.vcf))

for FILE in "${FILES[@]}"; do
# Extract the sample name from the filename
SAMPLE_NAME=$(basename "$FILE" .vcf)

# Define the output directory for the QC report
OUT_DIR_QC="$VCF_DIR/qc/$SAMPLE_NAME"

# Create the directory if it doesn't exist
mkdir -p $OUT_DIR_QC

# Calculate Ti/Tv ratio
bcftools stats $FILE > $OUT_DIR_QC/$SAMPLE_NAME.vcf.stats
plot-vcfstats -p $OUT_DIR_QC/plots/ $OUT_DIR_QC/$SAMPLE_NAME.vcf.stats

# Calculate missingness
java -jar $GATK -T VariantFiltration -R $REFERENCE -V $FILE --missingValuesInExpressionsShouldEvaluateAsFailing --filterExpression "vc.isNotFiltered()" --filterName "All Filters" -o $OUT_DIR_QC/$SAMPLE_NAME.filtered.vcf

# Calculate missingness per sample
java -jar $GATK -T SelectVariants -R $REFERENCE -V $OUT_DIR_QC/$SAMPLE_NAME.filtered.vcf -o $OUT_DIR_QC/$SAMPLE_NAME.selected.vcf

done
