#!/bin/bash

# 3) Convert to tab-delim, phased

# Columns: chr, position, hap1 allele, hap2 allele
# VCF if genotype info is 0/1 then unphased, if stored as 0|1 = phased (i.e hap1 | hap2)
# bcftools query -f '%CHROM\t%POS[\t%GT]\n' freeze10b.whi_only_chr21.bcf | head -3
# bcftools query f'[%CHROM\t%POS\t%SAMPLE\t%GT\n]  freeze10b.whi_only_chr21.bcf | head -3


# chr pos 0|1 >> need to find | in GT field, and replace delim with \t
DIR=/fh/scratch/delete90/kooperberg_c/mjohnson/cseqtl/results/genotype/split
OUTDIR=${DIR}/sub

for file in `cat geno_files_filt.txt`;
do cp $file $OUTDIR;
done

echo "File pasting complete"

cd $OUTDIR

for file in `cat ${DIR}/geno_files_filt.txt`;
do bcftools query -f '%CHROM\t%POS[\t%GT]\n' $file | sed 's:|:\t:g' -> ${file}_2;
done

echo "tab-delim haplotype files sorted"

# new header
