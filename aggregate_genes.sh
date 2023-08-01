#!/bin/bash

tree_dir=/fh/scratch/delete90/kooperberg_c/mjohnson/cseqtl/test_july/all_files/ASE
cd ${tree_dir}

# Get a list of all the *output.trecase.txt files in the current directory
files=( *output.trecase.txt )

# Create the header line for the output files
header="sample_id\ttotal\thap1\thap2\tASREC"

# Loop through each file
for file in "${files[@]}"; do
  # Extract the sample ID from the file name
  sample_id="${file%%-*}"

  # Debug output
  echo "Processing file: $file"
  echo "Sample ID: $sample_id"

  # Read the first ten rows of the file
  transcript_counts=$(awk -F'\t' 'NR > 1 && NR <= 11 { print $1 "\t" $2 "\t" $3 "\t"$4 "\t"$5 "\t" $3 + $4 }' "$file")

  # Create separate output files for each transcript
  for i in {1..10}; do
    # Extract the i-th transcript's counts
    transcript_line=$(echo "$transcript_counts" | awk -v i="$i" -F'\t' 'NR == i { print $0 }')

    # Extract the transcript ID and counts
    transcript=$(echo "$transcript_line" | cut -f1)
    trec=$(echo "$transcript_line" | cut -f2)
    hap1=$(echo "$transcript_line" | cut -f3)
    hap2=$(echo "$transcript_line" | cut -f4)
    hapN=$(echo "$transcript_line" | cut -f5)
    ASREC=$(echo "$transcript_line" | cut -f6)
 
    # Create the output file for the transcript
    output_file="${transcript}_counts.txt"   
    
    # Add header
    if [[ ! -f "$output_file" ]]; then
      echo -e "$header" > "$output_file"
    fi

    echo -e "$sample_id\t$trec\t$hap1\t$hap2\t$ASREC" >> "$output_file"
  done
done

# TODO add output dir or mv files