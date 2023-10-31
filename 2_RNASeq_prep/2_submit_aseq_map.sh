#!/bin/bash

#SBATCH --array=[3-1327]
#SBATCH --job-name=aseq_map
#SBATCH --nodes=1
#SBATCH --cpus-per-task=2
#SBATCH --mem-per-cpu=4G
#SBATCH --output=Rout/%j-%a.out
#SBATCH --error=Rerr/%j-%a.err
#SBATCH --mail-type=begin,end,fail
#SBATCH --mail-user=mjohnso5@fredhutch.org

ml fhR/4.2.2-foss-2021b

# Input
MANIFEST=lls_rna_array.csv
IN_DIR=/fh/scratch/delete90/kooperberg_c/mjohnson/cseqtl/results/genotype/LLS/Sept/LLS_Hap/aseq/ # Where aseq heterozygous snp files are stored
BAM_DIR=/fh/scratch/delete90/kooperberg_c/lls_rna/bam_files/bam_files/ # rnaseq bam files
REF_DIR=/fh/scratch/delete90/kooperberg_c/mjohnson/cseqtl/results/resources/ # Where exon_rds file is stored

# Output
OUT_DIR1=/fh/scratch/delete90/kooperberg_c/mjohnson/cseqtl/results/rnaseq/ # rnaseq mapped reads
OUT_DIR2=/fh/scratch/delete90/kooperberg_c/mjohnson/cseqtl/results/ASE/ # ASE mapped reads

echo "CPUs per task: $SLURM_CPUS_PER_TASK"
echo "Memory per CPU: $SLURM_MEM_PER_CPU"

Rscript 2_aseq_map.R ${MANIFEST} ${IN_DIR} ${BAM_DIR} ${REF_DIR} ${OUT_DIR1} ${OUT_DIR2} ${SLURM_ARRAY_TASK_ID}

# Original test subset = Ran in fh/scratch/delete90/kooperberg_c/mjohnson/cseqtl as sbatch -J manifest_rna_geno.csv 4_rnaseq_ase.sh
# LLS Run in fh/scratch/delete90/kooperberg_c/mjohnson/cseqtl/scripts/rnaseq_scripts as sbatch 2_submit_aseq_map.sh

# Manifest contains 3 cols: index, topmed_nwdid, lls_torid
# topmed_nwdid = wgs sample id (previously nwgc_sample_id)
# lls_torid = rna bam file prefix

