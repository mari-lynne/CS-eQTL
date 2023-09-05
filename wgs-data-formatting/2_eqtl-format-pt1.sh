#!/bin/bash

#SBATCH --mem=16
#SBATCH --output=Rout/%j-%a.out   # File to which STDOUT will be written, including job ID (%j)
#SBATCH --error=Rerr/%j-%a.err    # File to which STDERR will be written, including job ID (%j)

# 16 G x 23 CHR X 50 samples = 18400 G memory = 920 CPU

# This script aims to take bcf files (already split by CHR)
# Then split them by sample and make two files:
# One for ASeq counting (sample_hap.txt), and another with alleles recoded for CSeQTL models (sample_hap_map.txt).

# The final input for ASeq extractAsreads we need: (Per SAMPLE FILES)
# $1  $2  $3  $4  $5   $6
# CHR POS REF ALT HAP1 HAP2

# The intermediate input for CSeQTL_GS: (Per SAMPLE FILE)
# $1  $2  $3
# CHR POS SNP
# The final input for CSeQTL_GS we need: (Per SNP FILE)
# $1        $2
# Sample_ID SNP

# Load BCFtools module
ml BCFtools/1.14-GCC-11.2.0

# Read chromosome from the environment variable
chr=$chr

# Input directories and files
GENO_FILE=freeze10b.whi_only_chr
IN_DIR=/fh/scratch/delete90/kooperberg_c/mjohnson/cseqtl/results/genotype/LLS
# OUTPUT
OUT_DIR=${IN_DIR}/LLS_Hap

# Directory check:
  mkdir -p $OUT_DIR

ASSIGN=$(awk -F"," -v N=$SLURM_ARRAY_TASK_ID '
  NR==1 {for(i=1;i<=NF;i++) vars[i]=$i; next}
  $1 == N {for(i=1;i<=NF;i++) printf("%s=%s; ", vars[i], $i)}
' sample_array.csv)

# Make those assignments
eval $ASSIGN
SAMPLE_NAME=${topmed_nwdid}

# Execute bcftools to extract sample-specific genotype data and pipe it into awk
bcftools view -s ${SAMPLE_NAME} -m2 -M2 -v ${IN_DIR}/${GENO_FILE}${chr}.bcf \
  | bcftools query -f '%CHROM\t%POS\t%REF\t%ALT[\t%TGT\n]\n' \
  | awk -v outdir=${OUT_DIR} -v sample=${SAMPLE_NAME} -v chr=${chr} 'BEGIN {OFS="\t"} 
        # If it's the first row (Header), modify and print it
        NR==1 {
          $NF="HAP1\tHAP2";
          # Print modified header row to the new hap file
          print > (outdir "/" sample "_hap" ${chr} ".txt");
          # Print modified header row to the new map file (without HAP2)
          print $1, $2, $3, $4, $5 > (outdir "/" sample "_hap" ${chr} "_map.txt");
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
          print original "\t" a[2] > (outdir "/" sample "_hap" ${chr} ".txt");
          
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
          print sample, $1 "_" $2, $5 >> (${chr} "_temp_snp_data.txt");
        }'


