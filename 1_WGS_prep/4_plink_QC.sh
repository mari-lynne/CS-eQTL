#!/bin/bash

#SBATCH --nodes=1
#SBATCH --job-name=plink_format_final
#SBATCH --mem=20G
#SBATCH --cpus-per-task=12
#SBATCH --mail-type=begin,end,fail
#SBATCH --mail-user=mjohnso5@fredhutch.org
#SBATCH --output=plink_%j.out
#SBATCH --error=plink_%j.err 


# Filter/QC data using plink (HWE already set as well as DP/Q, MIND and GENO removed no variants)
# Reran with plink pfile to keep allele order 
# run script locally using sbatch in LLS_Hap dir for now, Plink2 in same folder

ml BCFtools
IN_FILE=whi_lls_concat.vcf.gz
OUT_FILE=whi_lls

# Convert bcf to plink ---------------------------------------------
./plink2 --vcf ${IN_FILE} --make-pgen --maf 0.01 --out ${OUT_FILE}

echo "VCF-PLINK conversion - DONE!"

# Convert final file back to BCF
./plink2 --pfile ${OUT_FILE} --export vcf --out ${OUT_FILE}_QC

# Add index
bcftools view ${OUT_FILE}_QC.vcf -Oz -o ${OUT_FILE}_QC.vcf.gz
bcftools index ${OUT_FILE}_QC.vcf.gz

echo "PLINK-VCF conversion - DONE!"

## Summary stats -----------------------------------------------------------

# TODO fix plot scripts

mv *.log log

