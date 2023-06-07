#!/bin/bash

ml BCFtools/1.14-GCC-11.2.0

### Run script in cseqtl/results as sbatch slurm_geno.sh

# 1) Concatenate chromosomes into bcf file
# 2) Split into sample (indv) specific bcf files
# 3) Convert to tab delim format for aseq

GENO_DIR=/fh/scratch/delete90/kooperberg_c/topmed_freeze10/phased/whi_only
OUT_DIR=/fh/scratch/delete90/kooperberg_c/mjohnson/cseqtl/results/genotype

INFILE=freeze10b.whi_only_chr
# OUTFILE=WHI

# OUT_DIR2=${OUT_DIR}/split
# OUT_DIR3=${OUT_DIR2}/sub

# 1) Concat --------------------------------------------------------------------

# Filter particpants across chromsomes
cd $OUT_DIR

for sample in ${GENO_DIR}${INFILE}{1..22}.bcf do;
bcftools view --samples-file ${OUT_DIR}/geno_files_filt.txt -Oz ${OUT_DIR}/${sample};
done

for chr in {1..22}; do
##         chr${chr}.P50.6147.unpruned.vcf.gz \
##         --samples-file ~/ratgtex/Eye/rat_ids.txt \
##         -Oz -o tmp_eye.chr${chr}.vcf.gz


# bcftools concat ${GENO_DIR}/freeze10b.whi_only_chr{1..22}.bcf --naive -Ob -o ${OUT_DIR}/${FILENAME}_concat.bcf
# 
# # 2) Split ---------------------------------------------------------------------
# 
# bcftools +split ${OUT_DIR}/${FILENAME}_concat.bcf -Ob -o ${OUT_DIR2}
# 
# # 3) Convert phased bcf file to tab-delim txt file for aseq
# # chr pos 0|1 >> need to find | in GT field, and replace delim with \t
# 
# for file in `cat geno_files_filt.txt`;
# do cp $file $OUTDIR;
# done
# 
# cd $OUTDIR
# 
# for file in `cat ${DIR}/geno_files_filt.txt`;
# do bcftools query -f '%CHROM\t%POS[\t%GT]\n' $file | sed 's:|:\t:g' -> ${file}_2;
# done