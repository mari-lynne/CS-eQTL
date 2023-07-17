#!/bin/bash

# genotype PC's
# Aims for SCT file set recaluclate genotype principle componenets using plink

# Slurm set up
#SBATCH --nodes=1
#SBATCH --mem-per-cpu=32000
#SBATCH --cpus-per-task=1
#SBATCH --output=plink_%j.out   # File to which STDOUT will be written, including job ID (%j)
#SBATCH --error=plink_%j.err    # File to which STDERR will be written, including job ID (%j)
#SBATCH --mail-type=begin
#SBATCH --mail-type=end
#SBATCH --mail-type=fail
#SBATCH --mail-user=mjohnso5@fredhutch.org

# Run script in genotype wd as sbatch wgs_pca.sh

ml plink/1.9-20200616

# Variables and Directories -----------------------------------------------------------

file=whi_concat_sub
dir=/fh/scratch/delete90/kooperberg_c/mjohnson/cseqtl/results/genotype/split
qcdir=${dir}/qc
pcadir=${qcdir}/pca

cd ${dir}

# Presuming these have already been pruned for IBD samples
# convert bcf to plink

plink \
--bcf ${file}.bcf \
--make-bed \
--out ${qcdir}/${file}

# Calculate LD

cd ${qcdir}

plink \
--bfile ${file} \
--indep-pairwise 50 5 0.2 \
--out ${file}.LD

# Prune variants
plink \
--bfile ${file} \
--allow-extra-chr \
--extract ${file}.LD.prune.in \
--make-bed \
--out ${file}.LD

# PCA on LD pruned, IBD and sib checked file
 
 plink \
 --bfile ${file}.LD \
 --pca \
 --out ${pcadir}/${file}.LD 

mv *.eigenval ${pcadir}
mv *.eigenvec ${pcadir}