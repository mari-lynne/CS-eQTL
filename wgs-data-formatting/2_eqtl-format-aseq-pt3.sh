#!/bin/bash

#SBATCH --mem=10000
#SBATCH --job-name=concat
#SBATCH --output=Rout/%j-%a.out   
#SBATCH --error=Rerr/%j-%a.err  
#SBATCH --mail-type=begin,end,fail
#SBATCH --mail-user=mjohnso5@fredhutch.org


# Directories
IN_DIR=/fh/scratch/delete90/kooperberg_c/mjohnson/cseqtl/results/genotype/LLS/LLS_Hap
OUT_DIR=${IN_DIR}/concat

# Create the output directory if it does not exist
mkdir -p $OUT_DIR

# Concatenate all hap files for this individual
cat ${IN_DIR}/${SAMPLE_NAME}_hap*.txt > ${OUT_DIR}/${SAMPLE_NAME}_all_hap.txt
