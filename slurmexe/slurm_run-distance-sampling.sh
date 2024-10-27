#!/bin/bash

# Check if a file name was provided
if [ -z "$1" ] || [ -z "$2" ]; then
  echo "Usage: $0 <cgraph-folder> <date>"
  exit 1
fi


cd ~/spanner/slurmexe
graphs=$1
date=$2
for file in $graphs/*; do
	echo $file;
	sbatch slurm_single-run-distance-sampling.slurm "$file" "$date"
done