#!/bin/bash
# RUN - run spanner algs (bash_convert-to-gbbs.sh)
cd ~/spanner/gbbs
snapfolder=~/spanner/graphs-sanitized-snap
gbbsfolder=~/spanner/graphs-sanitized-gbbs
for file in "$gbbsfolder"/*; do
	bazel run benchmarks/Spanner/FGV_Baseline:Spanner_main -- -s "$file"
	bazel run benchmarks/Spanner/MPVX_Baseline:Spanner_main -- -s "$file"
	bazel run benchmarks/Spanner/MPVX_CompactSpanner:Spanner_main -- -s "$file"
done