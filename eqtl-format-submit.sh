#!/bin/bash

# Initialize an empty array to hold the Slurm job IDs
job_ids=()

# Loop over each sample name and submit a job
while read -r SAMPLE_NAME; do
  # Submit the job and store the returned job ID in 'job_id'
  job_id=$(sbatch --parsable eqtl-format_pt1.sh ${SAMPLE_NAME})
  
  # Append the new job ID to the 'job_ids' array
  job_ids+=($job_id)
done < sample_list.txt

# Convert the job_ids array to a comma-separated string
job_ids_str=$(IFS=, ; echo "${job_ids[*]}")

# Submit the final job, specifying that it should only start after all the previous jobs have completed
final_job_id=$(sbatch --dependency=afterok:${job_ids_str} run_script_part2.sh)

# Optional: print the final job ID
echo "Final aggregation job submitted with ID: $final_job_id"
