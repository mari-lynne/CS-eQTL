#!/bin/bash

#SBATCH --job-name=ASeq_Format
#SBATCH --nodes=1
#SBATCH --mem-per-cpu=8G
#SBATCH --cpus-per-task=1
#SBATCH --output=Rout/%j-%a.out  
#SBATCH --error=Rerr/%j-%a.err 
#SBATCH --mail-type=begin,end,fail
#SBATCH --mail-user=mjohnso5@fredhutch.org
#SBATCH --array=[1-1327%500]

# Aims:
# Split reformated heterozygous bcf file into individual sample files
# Filter for heterozygous snps and reformat genotype column (|=tab delim)

# Info:
# Run in parallel jobs using slurm and sample array for efficiency
# Submit script as sbatch -J lls_sample_array_input.csv 3_aseq_format.sh

# 1) Variables and Directories ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

ml BCFtools/1.14-GCC-11.2.0

IN_DIR=/fh/scratch/delete90/kooperberg_c/mjohnson/cseqtl/results/genotype/LLS/Sept/LLS_Hap
IN_FILE=whi_lls_concat.bcf.vcf.gz
SAMPLE_ARRAY=lls_sample_array_input.csv

OUT_DIR=/fh/scratch/delete90/kooperberg_c/mjohnson/cseqtl/results/genotype/LLS/Sept/LLS_Hap/aseq

# Start time
start_time=$(date +"%Y-%m-%d %H:%M:%S")
echo "Task started at: $start_time"

# 2) Assign samples from csv array ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

SAMPLE_ASSIGN=$(awk -F"," -v N=$SLURM_ARRAY_TASK_ID '
  NR==1 {for(i=1;i<=NF;i++) vars[i]=$i; next}
  $1 == N {for(i=1;i<=NF;i++) printf("%s=%s; ", vars[i], $i)}
' $SAMPLE_ARRAY)

# Make assignments
eval $SAMPLE_ASSIGN
f=${topmed_nwdid}
echo "Sample ID = $f"

# 3) Make ASeq heterozygous SNP file ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# 1) Subset the bcf file by sample $(f)
# 2) Query relevant genotype fields [TGT]
# 3) Filter phased heterozygous genotype information (-i'GT="0|1" | GT="1|0"')
# 4) Reformat genotype field to tab delim and save output files

bcftools view -s ${f} ${IN_DIR}/${IN_FILE} \
      | bcftools query -f '%CHROM\t%POS[\t%TGT\n]\n' -i 'GT="0|1" | GT="1|0"' \
      | sed 's:|:\t:g' - > ${OUT_DIR}/${f}_hap.txt


echo "Hetrozygous SNP File - DONE!"

end_time=$(date +"%Y-%m-%d %H:%M:%S")
echo "Task ended at: $end_time"
