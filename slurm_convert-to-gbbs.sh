#!/bin/bash
snapfolder=~/spanner/graphs-sanitized-snap
for file in "$snapfolder"/*; do
	sbatch slurm_single-convert-to-gbbs.slurm "$file"
done