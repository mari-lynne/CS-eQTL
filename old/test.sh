#!/bin/bash
#
# Script to extract and modify WGS bcf files into aseq format

ml BCFtools/1.14-GCC-11.2.0

GENO_DIR=/fh/scratch/delete90/kooperberg_c/topmed_freeze10/phased/whi_only
OUT_DIR=/fh/scratch/delete90/kooperberg_c/mjohnson/cseqtl/results/genotype/split
FILE=freeze10b.whi_only_chr

# Steps ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# 1) Filter WHI particpants across chromsomes
# 2) Concatenate chromosomes into main bcf file
# 3) Split into sample (indv) specific bcf files
# 4) Convert to tab delim format for aseq

# 1) Filter particpants across chromsomes ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# 22 participants in geno_filt.txt, these have both RNAseq BAM files and WGS bcf files
# Use @ notation to subset 1 dim array

cd ${OUT_DIR}

chr=( 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X )
for f in "${chr[@]}"
do
  sbatch -c 1 -t 3-0 --mem 8000 --wrap="bcftools view -Ob --force-samples -S geno_files_filt.txt ${GENO_DIR}/${FILE}${f}.bcf > ${OUT_DIR}/${FILE}${f}.bcf"
done

# 2) Concatenate samples ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Run script using grabnode 20 GB mem, 1 CPU 

bcftools concat freeze10b.whi_only_chr{1..22}.bcf --naive -Ob -o whi_sub_concat.bcf

echo "concatenating DONE"

# for i in {1..22};do
#   bcftools index /fh/scratch/delete90/kooperberg_c/topmed_freeze10/phased/whi_only/freeze10b.whi_only_chr${i}.bcf
# done

# 3) Split indv samples ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Update later to do as sbatch array with samples to speed up probably

arr=( $(cat geno_files_filt.txt) )
echo ${arr}

for f in "${arr[@]}"
do
  sbatch -c 4 -t 2-0 --mem 20000 --mail-type=END --mail-type=FAIL --mail-user=mjohnso5@fredhutch.org --wrap="bcftools +split -S ${f} whi_sub_concat.bcf -Ob -o sub"
done

# 4) Split into tab format ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

cd ${OUT_DIR}/sub

for file in *.bcf; do
  bcftools query -f '%CHROM\t%POS[\t%GT]\n' ${file} | sed 's:|:\t:g' - > ${file}_2;
done


echo "tab-delim haplotype files sorted"




#!/bin/bash

arr=($(cat geno_files_filt.txt))
echo ${arr}
echo ${arr[3]}
echo "${arr[3]}"

arr2=$(cat geno_files_filt.txt)
echo "test 2"
echo ${arr2}
echo ${arr2[3]}
echo "${arr2[3]}"

arr3=`cat geno_files_filt.txt`
echo "test 3"
echo ${arr3}
echo ${arr3[3]}
echo "${arr3[3]}"



