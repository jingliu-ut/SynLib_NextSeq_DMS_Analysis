#!/bin/bash

# bash 06_01_dSortSeq_score_batch.sh -m <int> -M <int>

while getopts 'm:M:' opt; do
    case "$opt" in
        m)
          min="$OPTARG"
          ;;
        M)
          max="$OPTARG"
          ;;
    esac
done


for downsample_num in $(seq "$min" "$max")
do
    python3 06_01_dSortSeq_score.py -n "$downsample_num"
done




