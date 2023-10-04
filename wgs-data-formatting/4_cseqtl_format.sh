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
# BCF file, filtered in plink MAF <0.01, make_plink.sh script in LL_HAP dir
# TODO organise scripts

# 1) Set Directories and variables ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 

ml BCFtools/1.14-GCC-11.2.0
ml fhR

IN_DIR=/fh/scratch/delete90/kooperberg_c/mjohnson/cseqtl/results/genotype/LLS/Sept/LLS_Hap
REF_DIR=/fh/scratch/delete90/kooperberg_c/mjohnson/cseqtl/results/genotype/ref_data
OUT_DIR=/fh/scratch/delete90/kooperberg_c/mjohnson/cseqtl/results/genotype/LLS/Sept/gene_data

IN_FILE="lls_concat.bcf"
COORD_FILE="autosome_eqtl_${SLURM_ARRAY_TASK_ID}.txt"

# 2) Make CSeQTL SNP/Gene matricies ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Steps:
# For each line of the coord file i.e a gene;
# Extract the the chrom/start/end info - use that to subset bcf file
# Then recode genotype columns, save to a file named after the gene id
# Call R script to reshape data for CSeQTL input

while IFS=$'\t' read -r gene_id chrom start end ; do
bcftools view -r "chr${chrom}:${start}-${end}" ${IN_DIR}/${IN_FILE} \
| bcftools query -f '[%SAMPLE\t%ID\t%CHROM\t%POS\t%REF\t%ALT\t%TGT\n]' \
| awk -F "\t" -v OFS="\t" '
    {
        # Split the genotype field into an array based on "|"
        split($5, geno_array, "[|]");
        genotype_code = "5"; # Default/NA code
        
        # Determine the genotype code based on the alleles in the genotype field
        if (geno_array[1] == $3 && geno_array[2] == $3) {
            genotype_code = "0"; # Homozygous reference
        } else if (geno_array[1] == $3 && geno_array[2] == $4) {
            genotype_code = "1"; # Heterozygous
        } else if (geno_array[1] == $4 && geno_array[2] == $3) {
            genotype_code = "2"; # Heterozygous (order reversed)
        } else if (geno_array[1] == $4 && geno_array[2] == $4) {
            genotype_code = "3"; # Homozygous alternate
        }
        
        # Print the sample ID, position, and genotype code
        print $1, $2 ":" $3 ":" $4, genotype_code;
    }' >> ${OUT_DIR}/${gene_id}.txt
    Rscript 4_cseqtl_format.R ${gene_id} ${OUT_DIR}

done < "${REF_DIR}/${COORD_FILE}"

echo "Done with gene_snp files"
