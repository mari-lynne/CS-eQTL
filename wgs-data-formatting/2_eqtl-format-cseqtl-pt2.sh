#!/bin/bash

#SBATCH --job-name=format-pt2
#SBATCH --output=Rout/%j-%a.out   
#SBATCH --error=Rerr/%j-%a.err  
#SBATCH --mem=12G

# Ensure the chromosome variable is set
if [ -z ${chr} ]; then
    echo "Chromosome variable is not set. Exiting."
    exit 1
fi

# Ensure the input file exists
input_file="chr${chr}_temp_snp_data.txt"
if [ ! -f ${input_file} ]; then
    echo "Input file does not exist. Exiting."
    exit 1
fi

# Begin awk script to create per-SNP files.
# The SLURM script is already parallelized per individual.
# Using '>>' to append to each SNP file for cumulative data across individuals.

# Input Format: SAMPLE_ID CHR POS SNP (per Sample)
# Output Format: Sample_ID SNP (per SNP)

awk '{
  # Create the filename using chromosome and position identifiers.
  snp_file = "chr_" $2 "_" $3 ".txt"
  
  # Append SAMPLE ID and SNP genotype to the corresponding SNP file.
  print $1 "\t" $3 >> snp_file
}' ${input_file}
