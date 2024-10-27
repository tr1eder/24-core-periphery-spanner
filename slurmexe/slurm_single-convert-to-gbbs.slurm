#!/bin/bash
#SBATCH --job-name=convert-to-gbbs_%j.slurm                 # Job name
#SBATCH --output=slurmlog/convert-to-gbbs_%j.log            # Output file
#SBATCH --error=slurmerr/convert-to-gbbs_%j.err             # Error file
#SBATCH --ntasks=1                                          # Number of tasks (usually 1 for single commands)
#SBATCH --time=01:00:00                                     # Time limit (e.g., 1 hour)
#SBATCH --mem-per-cpu=64G                                   # Memory requirement
cd ~/spanner/gbbs
snapfolder=$1
gbbsfolder=$2
filename=$(basename $3)                                     # Get filename from SLURM_ARRAY_TASK_ID
tmpfile=$(mktemp)
~/spanner/convert_snap-bazel.sh "$snapfolder/$filename" > "$tmpfile"
bazel --output_base=~/spanner/slurmtmp/bazel_out_$SLURM_JOB_ID run //utils:snap_converter -- -s -i "$tmpfile" -o "$gbbsfolder/$filename"
rm "$tmpfile"