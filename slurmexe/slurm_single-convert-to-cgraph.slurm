#!/bin/bash
#SBATCH --job-name=convert-to-cgraph_%j.slurm       # Job name
#SBATCH --output=slurmlog/convert-to-cgraph_%j.log  # Output file
#SBATCH --error=slurmerr/convert-to-cgraph_%j.err   # Error file
#SBATCH --ntasks=1                                  # Number of tasks (usually 1 for single commands)
#SBATCH --time=01:00:00                             # Time limit (e.g., 1 hour)
#SBATCH --mem-per-cpu=128G                           # Memory requirement
cd ~/wormhole
file=$1
date=$2
filename=$(basename $file)
graphname="${filename%.edges}"
relpath=$(dirname "$(realpath --relative-to=$HOME/wormhole/graphs "$file")")
exe/createCGraph -g $graphname -i $relpath -o cgraphs-"$date"