#!/bin/bash

IN_DIR=/fh/scratch/delete90/kooperberg_c/mjohnson/cseqtl/results/genotype/LLS
SCRIPT_DIR=/fh/scratch/delete90/kooperberg_c/mjohnson/cseqtl/scripts
cd ${IN_DIR}

# Array of chromosomes
chromosomes=("1" "2")

# Loop through all chromosomes and submit jobs
for chr in "${chromosomes[@]}"; do
  # First, submit the inner sample-level jobs for this chromosome
  num_samples=$(($(wc -l < ${IN_DIR}/lls_sample_array_input.csv) - 1))  # Count the number of samples (subtract 1 for the header)
  
  # Submit the sample jobs and capture the main job ID
  first_batch_job=$(sbatch --array=1-${num_samples}%50 --export=chr=${chr} ${SCRIPT_DIR}/2_eqtl-format-pt1.sh | awk '{print $NF}')
  
  # Proceed aseq and cseqtl formatting scripts
  job2=$(sbatch \
  --parsable \
  --export=chr=$chr \
  --dependency=afterok:${first_batch_job} \
  ${SCRIPT_DIR}/2_eqtl-format-cseqtl-pt2.sh)

  job3=$(sbatch \
  --parsable \
  --array=1-${num_samples}%50 \
  --export=chr=$chr \
  --dependency=afterok:${first_batch_job} \
  ${SCRIPT_DIR}/2_eqtl-format-aseq-pt3.sh)

  # Print out the job IDs for logging/debugging
  echo "For chromosome $chr, submitted jobs dependent on sample-level jobs and their IDs are $job2, $job3"
done
