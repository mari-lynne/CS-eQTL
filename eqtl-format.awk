#!/bin/awk

# In prep make hap file to include ref/alt - we can exclude these fields later in R will be more efficient (so TODO update 2_wgs_data_prep)
# So lets presume we start with tab delim file 

# Dummy tab delim text file. Aim is to recode HAP1 and HAP2 alles to 0,1,2,3 based on if they match the ref or alt
# We need two files the first as it is for aseq counting, the second recoded to 0123 for the CSeQTL model. 
# Maybe see if we can combine this code with previous wgs code to do all in one.

# $1  $2  $3  $4  $5   $6
# CHR POS REF ALT HAP1 HAP2
# 1   123 A   T   A   A
# 1   243 C   G   G   C

# Output file 1 is
bcftools view -s ${SAMPLE_NAME} -m2 -M2 -v ${IN_DIR}/${GENO_FILE}${chr}.bcf \
  | bcftools query -f '%CHROM\t%POS\t%REF\t%ALT[\t%TGT\n]\n' \
  | awk 'BEGIN {OFS="\t"} 
        NR==1 {
          $NF="HAP1\tHAP2"; 
          print;
          next
        }
        {
          split($NF, a, "|");
          $NF=a[1]; 
          print $0 "\t" a[2]
        }' \
  > ${OUT_DIR2}/${SAMPLE_NAME}_hap${chr}.txt


# This makes this file

# $1  $2  $3  $4  $5   $6
# CHR POS REF ALT HAP1 HAP2
# 1   123 A   T   A   A
# 1   243 C   G   G   C


-F '\t' 'BEGIN {OFS = "\t"} \
NR == 1 {
  print $1, $2, $3, $4, $5;
  next
}
{
  if ($5 == $3 && $6 == $3) {
    $5 = "0";
  } else if ($5 == $3 && $6 == $4) {
    $5 = "1";
  } else if ($5 == $4 && $6 == $4) {
    $5 = "2";
  } else {
    $5 = "3";
  }
  print $1, $2, $3, $4, $5;
}' ${SAMPLE_NAME}_hap.txt > ${SAMPLE_NAME}_hap${chr}_map.txt


# Print modified header row (NR == 1)
# First if else checks hap1 and 2 fields for homozygous ref and recodes to 0 (AA)
# Second if else statement checks for het (hap2) = 1 (AB)
# Third if else checks for het (hap 1) = 2 (BA)
# Last else checks for homozygous alt  = 3 (BB)

# Combined version :)

bcftools view -s ${SAMPLE_NAME} -m2 -M2 -v ${IN_DIR}/${GENO_FILE}${chr}.bcf \
  | bcftools query -f '%CHROM\t%POS\t%REF\t%ALT[\t%TGT\n]\n' \
  | awk -v outdir=${OUT_DIR2} -v sample=${SAMPLE_NAME} -v chr=${chr} 'BEGIN {OFS="\t"} 
        NR==1 {
          $NF="HAP1\tHAP2";
          print > (outdir "/" sample "_hap" chr ".txt");
          print $1, $2, $3, $4, $5 > (outdir "/" sample "_hap" chr "_map.txt");
          next
        }
        {
          # Original HAP
          split($NF, a, "|");
          original=$0;
          $NF=a[1];
          print original "\t" a[2] > (outdir "/" sample "_hap" chr ".txt");
          
          # Recoded HAP
          if (a[1] == $3 && a[2] == $3) {
            $5 = "0";
          } else if (a[1] == $3 && a[2] == $4) {
            $5 = "1";
          } else if (a[1] == $4 && a[2] == $4) {
            $5 = "2";
          } else {
            $5 = "3";
          }
          print $1, $2, $3, $4, $5 > (outdir "/" sample "_hap" chr "_map.txt");
        }'
