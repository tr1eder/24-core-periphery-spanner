#!/bin/bash
#SBATCH --job-name=run-distance-sampling_%j.slurm       # Job name
#SBATCH --output=slurmlog/run-distance-sampling_%j.log  # Output file
#SBATCH --error=slurmerr/run-distance-sampling_%j.err   # Error file
#SBATCH --ntasks=1                                      # Number of tasks (usually 1 for single commands)
#SBATCH --time=01:00:00                                 # Time limit (e.g., 1 hour)
#SBATCH --mem-per-cpu=16G                               # Memory requirement
cd ~/wormhole
file=$1
date=$2
filename=$(basename $file)
csrgraphname="${filename%_CGraph.bin}" # eg. "graphname_spanner_FGVbase_CGraph.bin" -> "graphname_spanner_FGVbase"
graphname="${csrgraphname%_spanner_*}" # eg. "graphname_spanner_FGVbase" -> "graphname"
exe/sample -g $csrgraphname -gi cgraphs-"$date" -i inputs-"$date"/"$graphname"_input.txt -o dists-"$date"
