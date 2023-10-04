#!/bin/bash

#SBATCH --mem=16G
#SBATCH --job-name
#SBATCH --output=Rout/%j-%a.out   
#SBATCH --error=Rerr/%j-%a.err

# 16 G x 23 CHR X 50 samples = 18400 G memory = 920 CPU

# This script aims to take bcf files (already split by CHR in the submit script) and split them by sample
# For each individual we need to make two files:
# One for ASeq counting (sample_hap.txt), and another with alleles recoded for CSeQTL models (sample_hap_map.txt).

# The final input for ASeq counting - extractAsreads: (Per SAMPLE FILES)
# $1  $2    $3   $4
# CHR POS  HAP1 HAP2

# The intermediate input for CSeQTL_GS: (Per SAMPLE FILE)
# $1  $2  $3
# CHR POS SNP

# The final input for CSeQTL_GS we need: (Per SNP FILE)
# $1        $2
# Sample_ID SNP

# With SNP recoded based on the ref/alt combination:
# "AA","AB","BA","BB" : 0,1,2,3

# Function to check if the previous command was successful
trap 'last_command=$BASH_COMMAND' DEBUG

check_command_success() {
    if [ $? -ne 0 ]; then
        echo "Command failed: $last_command"
        exit 1
    fi
}

# Retrieve chromosome from environment variable
chr=${chr:-"NotSet"}
if [ "$chr" == "NotSet" ]; then
    echo "Chromosome not set. Exiting."
    exit 1
fi

# Input directories and files
GENO_FILE=freeze10b.whi_only_chr
IN_DIR=/fh/scratch/delete90/kooperberg_c/mjohnson/cseqtl/results/genotype/LLS
# OUTPUT
OUT_DIR=${IN_DIR}/LLS_Hap

# Extract sample ids from slurm array
SAMPLE_ASSIGN=$(awk -F "," -v N=$SLURM_ARRAY_TASK_ID '
  NR==1 {for(i=1;i<=NF;i++) vars[i]=$i; next}
  $1 == N {for(i=1;i<=NF;i++) printf("%s=%s; ", vars[i], $i)}
' lls_sample_array_input.csv)
check_command_success

# Make those assignments
eval $SAMPLE_ASSIGN
SAMPLE_NAME=${topmed_nwdid}

# Temp file
TEMP_SNP_DATA="${chr}_temp_snp_data.txt"

echo "$SAMPLE_NAME"

# Load BCFtools module
ml BCFtools/1.14-GCC-11.2.0

# Main code - make aseq and cseqtl intermediate files from chr bcf files ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# -m2 -M2 -v snps selects only biallelic snps
bcftools view -s ${SAMPLE_NAME} -m2 -M2 -v snps ${IN_DIR}/${GENO_FILE}${chr}.bcf |
bcftools query -f '%CHROM\t%POS\t%REF\t%ALT[\t%TGT\n]\n' |
awk -v outdir=${OUT_DIR} -v sample=${SAMPLE_NAME} -v chr=${chr} -v temp_snp_data=${TEMP_SNP_DATA} 'BEGIN {OFS="\t"} \
{
  # 1) Begin recoding bcf files for ASeq ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  # Split the last field (genotype) into array 'geno_array' separated by the phase delimiter "|"
  split($NF, geno_array, "|");

  # Print modified line: CHR POS HAP1 HAP2
  print $1, $2, geno_array[1], geno_array[2] > (outdir "/" sample "_hap_" chr ".txt");

  # 2) Begin recoding for CSeQTL model ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  # Skip if any crucial field is missing
  if ($1 == "" || $2 == "" || $3 == "" || $4 == "" || $5 == "") {
    next;
  }
  
  # Initialize variable to hold the genotype code (0, 1, 2, or 3)
  genotype_code = "";

  # Determine genotype code based on the ref/alt condition
  if (geno_array[1] == $3 && geno_array[2] == $3) {
    genotype_code = "0";
  } else if (geno_array[1] == $3 && geno_array[2] == $4) {
    genotype_code = "1";
  } else if (geno_array[1] == $4 && geno_array[2] == $3) {
    genotype_code = "2";
  } else if (geno_array[1] == $4 && geno_array[2] == $4) {
    genotype_code = "3";
  } else {
    genotype_code = "NA";
  }

  # Print the output as SAMPLE_ID, CHR, POS, GENOTYPE_CODE
  print sample, $1, $2, genotype_code >> (outdir "/chr" temp_snp_data)
}'
check_command_success

# Next steps:
# Concatenate aseq files by chromsome
# Concatenate cseqtl files by SNP
