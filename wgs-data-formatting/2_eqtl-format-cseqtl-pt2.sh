#!/bin/bash

#SBATCH --job-name=format-pt2
#SBATCH --output=Rout/%j-%a.out   
#SBATCH --error=Rerr/%j-%a.err  
#SBATCH --mem=12G

# Ensure the chromosome variable is set
if [ -z "$chr" ]; then
    echo "Chromosome variable is not set. Exiting."
    exit 1
fi

# Ensure the input file exists
input_file="${chr}_temp_snp_data.txt"
if [ ! -f "$input_file" ]; then
    echo "Input file does not exist. Exiting."
    exit 1
fi

awk '{
  snp_file = $2 ".txt"
  print $1 "\t" $3 >> snp_file
}' ${chr}_temp_snp_data.txt