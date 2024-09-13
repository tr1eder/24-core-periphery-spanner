#!/bin/bash
gbbsfolder=~/spanner/graphs-sanitized-gbbs
for file in "$gbbsfolder"/*; do
	sbatch slurm_single-run-spanners.slurm "$file" "FGV_Baseline"
	sbatch slurm_single-run-spanners.slurm "$file" "MPVX_Baseline"
	sbatch slurm_single-run-spanners.slurm "$file" "MPVX_CompactSpanner"
done