#!/bin/bash
# gbbsfolder=~/spanner/graphs-sanitized-gbbs$

# Check if the usage is correct
if [ -z "$1" ]; then
	echo "Usage: $0 <gbbs-folder>"
	exit 1
fi

gbbsfolder=$1
for file in $gbbsfolder/*; do
	echo $file;
	sbatch slurm_single-run-spanners.slurm "$file" "FGV_Baseline"
	sbatch slurm_single-run-spanners.slurm "$file" "MPVX_Baseline"
	sbatch slurm_single-run-spanners.slurm "$file" "MPVX_CompactSpanner"
done