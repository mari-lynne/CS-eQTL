#!/bin/bash
set -e

#SBATCH --nodes=1
#SBATCH --mem-per-cpu=20G
#SBATCH --cpus-per-task=6
#SBATCH --output=pca_%j.out
#SBATCH --error=pca_%j.err
#SBATCH --mail-type=begin,end,fail
#SBATCH --mail-user=mjohnso5@fredhutch.org

# Variables and Directories -----------------------------------------------------------

file=whi_lls

# Remove duplicate variants
./plink2 \
--pfile ${file} \
--set-missing-var-ids @:#\$r:\$a \
--snps-only just-acgt \
--rm-dup exclude-all \
--make-pgen \
--out ${file}.QC

# 1) Prune data ------------------------------------------------------------------------

./plink2 \
--pfile ${file}.QC \
--indep-pairwise 50 5 0.2 \
--out ${file}.LD

# Prune variants
./plink2 \
--pfile ${file}.QC \
--extract ${file}.LD.prune.in \
--make-pgen \
--out ${file}.LD

# 2) Calculate IBD --------------------------------------------------------------------

./plink2 \
--pfile ${file}.LD \
--make-king \
--out ${file}.LD

# Remove one sample from each pair with pi-hat (% IBD) above threshold (0.1875):
#awk '$10 >= 0.1875 {print $2, $4, $10}' ${file}.LD.genome > pair_list.txt
#awk '$10 >= 0.1875 {print $1, $2}' ${file}.LD.genome | uniq > ${file}.outliers.txt
#wc -l ${file}.outliers.txt

#echo "Outlier list - Done"

# Filter IBD outliers in PLINK

#./plink2 \
#--pfile ${file}.LD \
#--remove ${file}.outliers.txt \
#--make-pgen \
#--out ${file}.IBD

# 3) PCA ------------------------------------------------------------------------------

 ./plink2 \
 --pfile ${file}.LD \
 --pca \
 --out ${file}.LD

