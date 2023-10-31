#!/bin/bash
#SBATCH --nodes=1
#SBATCH --job-name=cseqtl_format
#SBATCH --mem-per-cpu=20G
#SBATCH --output=Rout/%j-%a.out  
#SBATCH --error=Rerr/%j-%a.err    
#SBATCH --mail-type=begin,end,fail
#SBATCH --mail-user=mjohnso5@fredhutch.org
#SBATCH --array=1-22

# Aims:
# Make recoded SNP genotype files per gene for CSeQTL

# Info:
# Array runs 22 instances of job, also used to select coordinate file
# Coordinate file made in R (Local/git dir/gene_coord_format.R) and subsetted by chromosome
# BCF file, filtered in 4_plink_QC MAF <0.01, data saved in LL_HAP dir

# 1) Set Directories and variables ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 

ml BCFtools/1.18-GCC-12.2.0
ml fhR

IN_DIR="/fh/scratch/delete90/kooperberg_c/mjohnson/cseqtl/results/genotype/LLS/Sept/LLS_Hap"
REF_DIR="/fh/scratch/delete90/kooperberg_c/mjohnson/cseqtl/results/genotype/ref_data"
OUT_DIR="/fh/scratch/delete90/kooperberg_c/mjohnson/cseqtl/results/genotype/LLS/Sept/gene_data/chr${SLURM_ARRAY_TASK_ID}"

IN_FILE="whi_lls_QC.vcf.gz"
COORD_FILE="${REF_DIR}/autosome_eqtl_${SLURM_ARRAY_TASK_ID}.txt"

# Checks:

if [ ! -d "$REF_DIR" ] || [ ! -d "$IN_DIR" ] || [ ! -d "$OUT_DIR" ]; then
  echo "One or more directories do not exist. Exiting."
  exit 1
fi

if [ ! -f "$COORD_FILE" ] || [ ! -f "${IN_DIR}/${IN_FILE}" ]; then
  echo "One or more files do not exist. Exiting."
  exit 1
fi

# 2) Make CSeQTL SNP/Gene matricies ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Steps:
# For each line of the coord file i.e a gene;
# Extract the the chrom/start/end info - use that to subset bcf file
# Then recode genotype columns, save to a file named after the gene id
# Call R script to reshape data for CSeQTL input

# 2) Make CSeQTL SNP/Gene matricies ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# NOTE! check chromsome annotation, if it's chr21 in vcf file then run view as chr${chrom}

while IFS=$'\t' read -r gene_id chrom start end ; do
bcftools view -r "${chrom}:${start}-${end}" ${IN_DIR}/${IN_FILE} |
bcftools query -f '[%SAMPLE\t%ID\t%CHROM\t%POS\t%REF\t%ALT\t%TGT\n]'|
awk 'BEGIN {FS="\t"; OFS="\t"} {split($7, geno_array, "[|]");
        genotype_code = "5"; 
        if (geno_array[1] == $5 && geno_array[2] == $5) {
            genotype_code = "0"; 
        } else if (geno_array[1] == $5 && geno_array[2] == $6) {
            genotype_code = "1";
        } else if (geno_array[1] == $6 && geno_array[2] == $5) {
            genotype_code = "2"; 
        } else if (geno_array[1] == $6 && geno_array[2] == $6) {
            genotype_code = "3"; 
        }
        print $1, $3 ":" $4 ":" $2, genotype_code;
    }' > ${OUT_DIR}/${gene_id}_temp.txt
Rscript 4_cseqtl_format.R ${gene_id} ${OUT_DIR}
done < ${COORD_FILE}

echo "CSeQTL SNP per Gene files - DONE!"

# Genotype codes:
# Hom ref = 0
# Het ref/alt = 1
# Het alt/ref = 2
# Hom alt = 3
# NA = 5
