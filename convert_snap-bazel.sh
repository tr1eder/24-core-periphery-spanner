#!/bin/bash

# Check if a file name was provided
if [ -z "$1" ]; then
  echo "Usage: $0 <graph-file>"
  exit 1
fi

graph_file="$1"

# Output the first line with a '#' prepended, followed by the rest of the file
{
  echo "# $(head -n 1 "$graph_file")"
  tail -n +2 "$graph_file"
}
