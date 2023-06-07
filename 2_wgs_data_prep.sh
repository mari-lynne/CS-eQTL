#!/bin/bash

#SBATCH --array=[1-20]
#SBATCH --nodes=1
#SBATCH --mem-per-cpu=32000
#SBATCH --cpus-per-task=1
#SBATCH --output=Rout/par-%j-%a.out   # File to which STDOUT will be written, including job ID (%j)
#SBATCH --error=Rout/par-%j-%a.out    # File to which STDERR will be written, including job ID (%j)
#SBATCH --mail-type=begin
#SBATCH --mail-type=end
#SBATCH --mail-type=fail
#SBATCH --mail-user=mjohnso5@fredhutch.org

# Script to extract and modify WGS bcf files into aseq format
# Run script in genotype wd as sbatch --test-only -J samples.csv split_indv.sh

# Steps ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# 1) Filter WHI particpants across chromsomes
# 2) Concatenate chromosomes into main bcf file
# 3a) Split into sample (indv) specific bcf files
# 3b) Convert to tab delim format for aseq


ml BCFtools/1.14-GCC-11.2.0

# Assign variables ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 

GENO_DIR=/fh/scratch/delete90/kooperberg_c/topmed_freeze10/phased/whi_only
OUT_DIR=/fh/scratch/delete90/kooperberg_c/mjohnson/cseqtl/results/genotype/split
FILE=freeze10b.whi_only_chr
OUT_DIR2=${OUT_DIR}/sub


ASSIGN=$(awk -F"," -v N=$SLURM_ARRAY_TASK_ID '
  NR==1 {for(i=1;i<=NF;i++) vars[i]=$i; next}
  $1 == N {for(i=1;i<=NF;i++) printf("%s=%s; ", vars[i], $i)}
' sample_array.csv)

# make those assignments
eval $ASSIGN

f=${file_prefix}

# 1) Filter particpants across chromsomes ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# 21 participants in geno_filt.txt, these have both RNAseq BAM files and WGS bcf files
# Use @ notation to subset 1 dim array

cd ${OUT_DIR}

chr=( 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X )
for f in "${chr[@]}"
do
  sbatch -c 1 -t 3-0 --mem 8000 --wrap="bcftools view -Ob --force-samples -S geno_files_filt.txt ${GENO_DIR}/${FILE}${f}.bcf > ${OUT_DIR}/${FILE}${f}.bcf"
done

echo "Sample Filtering DONE"

# 2) Concatenate samples ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

bcftools concat freeze10b.whi_only_chr{1..22}.bcf --naive -Ob -o whi_sub_concat.bcf

echo "Chr Concatenating DONE"

# 3) Split bcf file by subset and extract phased genotype information ~~~~~~~~~~

bcftools view -c1 -s ${f} whi_sub_concat.bcf \
        | bcftools query -f '%CHROM\t%POS[\t%TGT\n]\n' \
        | sed 's:|:\t:g' - > ${OUT_DIR2}/${f}_hap.txt

# BCFtools options:
# -c1 = keep only variants with at least one alternate allele (as is necessary for ASE)
# -s ${f} select sample in bcf concat file from array var
# -Ob -o = output bcf files
#  query -f select fields from bcf file, tab delim
#  sed replace 0|1 genotype field with tab delim for aseq2

# bcftools view -c1 -s ${f} whi_sub_concat.bcf -Ob -o ${f}.bcf

echo "File Splitting complete"
