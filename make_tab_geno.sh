#!/bin/bash

# 3) Convert phased bcf file to tab-delim phased file for aseq:
# Columns: chr, position, hap1 allele, hap2 allele
# Use bcf tools query to output fields, and sed to find and replace delim

# bcftools query -f '%CHROM\t%POS[\t%GT]\n' freeze10b.whi_only_chr21.bcf | head -3
# chr pos 0|1 >> need to find | in GT field, and replace delim with \t

ml BCFtools/1.14-GCC-11.2.0

DIR=/fh/scratch/delete90/kooperberg_c/mjohnson/cseqtl/results/genotype/split
OUTDIR=${DIR}/sub

for file in `cat geno_files_filt.txt`;
do cp $file $OUTDIR;
done

echo "File pasting complete"

cd $OUTDIR

for file in `cat ${DIR}/geno_files_filt.txt`;
do bcftools query -f '%CHROM\t%POS[\t%GT]\n' $file | sed 's:|:\t:g' -> ${file};
done

echo "tab-delim haplotype files sorted"
