#!/bin/bash

# Usage: bash 02_downsample.sh -m <min> -M <max>
# -m: first number of downsampled files
# -M: last number of downsampled files


# Parse command-line arguments for m and M
while getopts 'm:M:' opt; do
    case "$opt" in
        m)
          min="$OPTARG" ;;
        M)
          max="$OPTARG" ;;
    esac
done

# Check if both min and max are set
if [ -z "$min" ] || [ -z "$max" ]; then
    echo "Error: -min and -max options are required."
    exit 1
fi

# Generate a for loop for min to max
for i in $(seq "$min" "$max")
do
  # Run the 02_downsample.py script with the random_state i
  echo "Running 02_downsample.py with random_state $i"
  python3 02_downsample.py --random_state "$i"
done

