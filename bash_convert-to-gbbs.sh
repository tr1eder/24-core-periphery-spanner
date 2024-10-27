#!/bin/bash
# RUN - convert from snap to gbbs format
cd ~/spanner/gbbs
snapfolder=~/spanner/graphs-sanitized-snap
gbbsfolder=~/spanner/graphs-sanitized-gbbs
for file in "$snapfolder"/*; do     
	filename=$(basename "$file");  
	# sed -i '1s/^/# /' "$file";
	bazel run //utils:snap_converter -- -s -i <(~/spanner/convert_snap-bazel.sh "$file") -o "$gbbsfolder/$filename"; 
done    