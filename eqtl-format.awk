#!/bin/bash

# This script aims to recode HAP1 and HAP2 alleles to 0,1,2,3 based on whether they match the ref or alt alleles.
# Two files are created: one for ASeq counting (sample_hap.txt) and another recoded for CSeQTL models (sample_hap_map.txt).

# The input for ASeq extractAsreads: (Per SAMPLE FILES)
# $1  $2  $3  $4  $5   $6
# CHR POS REF ALT HAP1 HAP2

# The intermediate input for CSeQTL_GS: (Per SAMPLE FILE)
# $1  $2  $3
# CHR POS SNP
# The final input for CSeQTL_GS: (Per SNP FILE)
# $1        $2
# Sample_ID SNP

# Execute bcftools to extract sample-specific genotype data and pipe it into awk
bcftools view -s ${SAMPLE_NAME} -m2 -M2 -v ${IN_DIR}/${GENO_FILE}${chr}.bcf \
  | bcftools query -f '%CHROM\t%POS\t%REF\t%ALT[\t%TGT\n]\n' \
  | awk -v outdir=${OUT_DIR2} -v sample=${SAMPLE_NAME} -v chr=${chr} 'BEGIN {OFS="\t"} 
        # If it's the first row (Header), modify and print it
        NR==1 {
          $NF="HAP1\tHAP2";
          # Print modified header row to the new hap file
          print > (outdir "/" sample "_hap" chr ".txt");
          # Print modified header row to the new map file (without HAP2)
          print $1, $2, $3, $4, $5 > (outdir "/" sample "_hap" chr "_map.txt");
          # Skip to the next row
          next
        }
        {
          # Split the last field (genotype) into array 'a' separated by "|"
          split($NF, a, "|");
          # Keep a copy of the original line
          original=$0;
          # Replace the last field (genotype) with HAP1 value
          $NF=a[1];
          # Print original line with HAP1 and new HAP2
          print original "\t" a[2] > (outdir "/" sample "_hap" chr ".txt");
          
          # Begin recoding for CSeQTL model ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
          # Homozygous ref allele
          if (a[1] == $3 && a[2] == $3) {
            $5 = "0";
          } 
          # Heterozygous, with ref allele first
          else if (a[1] == $3 && a[2] == $4) {
            $5 = "1";
          } 
          # Heterozygous, with alt allele first 
          else if (a[1] == $4 && a[2] == $3) {
            $5 = "2";
          } 
          # Homozygous alt allele (or any other combination not accounted for)
          else {
            $5 = "3";
          }
          # Print the recoded data for the CSeQTL model
          print sample, $1 "_" $2, $5 >> (chr "temp_snp_data.txt");
;
        }'

# The second phase of the script needs to aggregate data per SNP accross individuals:

awk '{
  snp_file = $2 ".txt"
  print $1 "\t" $3 >> snp_file
}' temp_snp_data.txt
