#!/bin/bash

# cd to trecase dir, save output to ciber dir
trecase_dir=~/Documents/CSeQTL/data/ciber_ase

cd ${trecase_dir}

# 1) get total counts from all files - save as a temporary version
for f in *trecase.txt; do
  awk '{print $2}' ${f} > ${f}_temp
done

# 2) paste total counts across all the samples
paste *_temp > temp_counts.txt

# 3a) Get 1st column from orginal file
input_file=`find *output.trecase.txt -type f -printf "%T@ %p\n" | sort -n | cut -d' ' -f 2- | head -n 1`

# 3b) paste gene column to the toal counts file to make final version 
cut $input_file -f 1 | paste - temp_counts.txt > total_counts.txt

# 4) remove temp files
rm *temp*
  