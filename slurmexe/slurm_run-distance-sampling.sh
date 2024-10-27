#!/bin/bash

# Check if a file name was provided
if [ -z "$1" ] || [ -z "$2" ]; then
  echo "Usage: $0 <cgraph-folder> <date>"
  exit 1
fi


cd ~/spanner/slurmexe
graphs=$1
date=$2

# substrings=("twitter2" "hollywood" "sinaweibo")
substrings=("soc-twitter-og")

for file in $graphs/*; do
    # Check if the filename contains any of the specified substrings
    for substring in "${substrings[@]}"; do
        if [[ "$file" == *"$substring"* ]]; then
            echo "$file"
            sbatch slurm_single-run-distance-sampling.slurm "$file" "$date"
            break  # Exit the inner loop once a match is found
        fi
    done
done