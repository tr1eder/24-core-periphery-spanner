#!/bin/bash

# Check if a file name was provided
if [ -z "$1" ]; then
  echo "Usage: $0 <graph-file>"
  exit 1
fi

graph_file="$1"

# Extract node and edge counts from the second line
node_count=$(sed -n '2p' "$graph_file" | awk '{print $3}')
edge_count=$(sed -n '2p' "$graph_file" | awk '{print $5}')

# Output node and edge counts followed by the rest of the file
{
  echo "$node_count $edge_count"
  tail -n +4 "$graph_file"
}
