#!/bin/bash

#SBATCH --nodes=1
#SBATCH --job-name=LLS_concat
#SBATCH --mem-per-cpu=20G
#SBATCH --cpus-per-task=12
#SBATCH --output=Rout/%j-%a.out  
#SBATCH --error=Rerr/%j-%a.err    
#SBATCH --mail-type=begin,end,fail
#SBATCH --mail-user=mjohnso5@fredhutch.org

# Aims:
# Concatenate chromosomes and make heterozygous SNP file (multisample) for aseq

ml BCFtools/1.14-GCC-11.2.0

# 1) Set Directories and variables ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 

IN_DIR=/fh/scratch/delete90/kooperberg_c/mjohnson/cseqtl/results/genotype/LLS/Sept
OUT_DIR=/fh/scratch/delete90/kooperberg_c/mjohnson/cseqtl/results/genotype/LLS/Sept/LLS_Hap

IN_FILE=freeze10b.whi_only_chr
OUT_FILE=whi_lls_reconcat
# reran 0ct 9th without GT filter

# 2) Concatenate Chromsomes ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

bcftools concat ${IN_DIR}/${IN_FILE}{1..22}.bcf -Oz -o ${OUT_DIR}/${OUT_FILE}.vcf.gz
echo "Chr Concatenating - DONE!"

bcftools index ${OUT_DIR}/${OUT_FILE}.vcf.gz
echo "File indexing - DONE!"

