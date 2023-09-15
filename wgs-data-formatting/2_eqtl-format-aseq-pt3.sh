#!/bin/bash

#SBATCH --mem=12G
#SBATCH --output=Rout/%j-%a.out   
#SBATCH --error=Rerr/%j-%a.err  
#SBATCH --mail-type=begin,end,fail
#SBATCH --mail-user=mjohnso5@fredhutch.org

# Aims concatenate haplotype files which have been split per chr/individual, into files grouped per ind (with multiple chromsomes)

# Extract sample assignments from array
SAMPLE_ASSIGNMENT=$(awk -F"," -v N=$SLURM_ARRAY_TASK_ID '
  NR==1 {for(i=1;i<=NF;i++) vars[i]=$i; next}
  $1 == N {for(i=1;i<=NF;i++) printf("%s=%s; ", vars[i], $i)}
' lls_sample_array_input.csv)

# Make those assignments
eval ${SAMPLE_ASSIGN}
SAMPLE_NAME=${topmed_nwdid}

# Check if SAMPLE_NAME is set
if [ -z ${SAMPLE_NAME} ]; then
    echo "SAMPLE_NAME is not set. Exiting."
    exit 1
fi

# Directories
IN_DIR=/fh/scratch/delete90/kooperberg_c/mjohnson/cseqtl/results/genotype/LLS/LLS_Hap
OUT_DIR=${IN_DIR}/concat

# Create the output directory if it does not exist
mkdir -p $OUT_DIR

# Check if any hap files exist for this individual
HAP_FILES=$(ls ${IN_DIR}/${SAMPLE_NAME}_hap*.txt 2> /dev/null)
if [ -z ${HAP_FILES} ]; then
    echo "No haplotype files found for individual ${SAMPLE_NAME}. Exiting."
    exit 1
fi

# Concatenate all hap files for this individual
cat ${HAP_FILES} > ${OUT_DIR}/${SAMPLE_NAME}_all_hap.txt


