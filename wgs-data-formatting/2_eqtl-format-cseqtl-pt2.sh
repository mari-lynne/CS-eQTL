#!/bin/bash

#SBATCH --job-name=format-pt2
#SBATCH --mem=16G

# Read chromosome from the environment variable
chr=$chr

awk '{
  snp_file = $2 ".txt"
  print $1 "\t" $3 >> snp_file
}' ${chr}_temp_snp_data.txt