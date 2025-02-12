#!/bin/bash
# Parse input arguments for gbbsfolder, filespecifier, spanners, spannertranslations
gbbsfolder=~/spanner/graphs-sanitized-gbbs
filespecifier=""
spanners=("fgvb" "fgvcs" "mpvxb" "mpvxcs")
declare -A spannertranslations=(
	["fgvb"]="FGV_Baseline"
	["fgvcs"]="FGV_CompactSpanner"
	["mpvxb"]="MPVX_Baseline"
	["mpvxcs"]="MPVX_CompactSpanner"
)

while [ $# -gt 0 ]; do
    case "$1" in
        -gbbsfolder)
            gbbsfolder="$2"
            shift 2
            ;;
        -filespecifier)
            filespecifier="$2"
            shift 2
            ;;
        -spanners)
            IFS=',' read -ra spanners <<< "$2"
            shift 2
            ;;
        *)
            echo "Usage: $0 [-gbbsfolder <folder>] [-filespecifier <specifier>] [-spanners <fgvb,fgvcs,mpvxb,mpvxcs,...>]"
            exit 1
            ;;
    esac
done

# Change directory into the gbbs folder as specified
# cd "$gbbsfolder"
cd ~/spanner/gbbs

for file in "$gbbsfolder/$filespecifier"*; do
	for i in "${!spanners[@]}"; do
		echo "===> Running ${spanners[i]} on $file"
		bazel run benchmarks/Spanner/${spannertranslations[${spanners[i]}]}:Spanner_main -- -s "$file"
	done
done

# bazel run benchmarks/Spanner/FGV_Baseline:Spanner_main -- -s "$file"
# bazel run benchmarks/Spanner/MPVX_Baseline:Spanner_main -- -s "$file"
# bazel run benchmarks/Spanner/MPVX_CompactSpanner:Spanner_main -- -s "$file"