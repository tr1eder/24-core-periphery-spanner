#!/bin/bash
#SBATCH --job-name=run-spanners_%j.slurm       # Job name
#SBATCH --output=slurmlog/run-spanners_%j.log  # Output file
#SBATCH --error=slurmerr/run-spanners_%j.err   # Error file
#SBATCH --ntasks=1                             # Number of tasks (usually 1 for single commands)
#SBATCH --time=01:00:00                        # Time limit (e.g., 1 hour)
#SBATCH --mem-per-cpu=8G                       # Memory requirement
cd ~/spanner/gbbs
file=$1
spannerversion=$2
bazel --output_base=~/spanner/slurmtmp/bazel_out_$SLURM_JOB_ID run benchmarks/Spanner/$spannerversion:Spanner_main -- -s "$file"