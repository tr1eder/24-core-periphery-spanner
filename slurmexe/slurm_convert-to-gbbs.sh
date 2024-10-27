#!/bin/bash

# Check if the usage is correct
if [ -z "$1" ] || [ -z "$2" ]; then
	echo "Usage: $0 <snap-folder> <gbbs-folder>"
	exit 1
fi

snapfolder=$1
gbbsfolder=$2
for file in $snapfolder/*; do
	echo $file;
	sbatch slurm_single-convert-to-gbbs.slurm "$snapfolder" "$gbbsfolder" "$file"
done