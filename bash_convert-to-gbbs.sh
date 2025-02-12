#!/bin/bash
# RUN - convert from snap to gbbs format

snapfolder=~/spanner/graphs-sanitized-snap
gbbsfolder=~/spanner/graphs-sanitized-gbbs
filespecifier=""

while [ $# -gt 0 ]; do
    case "$1" in
        -snapfolder)
            snapfolder="$2"
            shift 2
            ;;
        -gbbsfolder)
            gbbsfolder="$2"
            shift 2
            ;;
        -filespecifier)
            filespecifier="$2"
            shift 2
            ;;
        *)
            echo "Usage: $0 [-snapfolder <snapfolder>] [-gbbsfolder <gbbsfolder>] [-filespecifier <filespecifier>]"
            exit 1
            ;;
    esac
done

cd ~/spanner/gbbs

for file in "$snapfolder/$filespecifier"*; do
    filename=$(basename "$file")
    tmpfile=$(mktemp)
    ~/spanner/convert_snap-bazel.sh "$file" > $tmpfile
    echo "===> Converting $filename"
    bazel run //utils:snap_converter -- -s -i $tmpfile -o "$gbbsfolder/$filename"
done