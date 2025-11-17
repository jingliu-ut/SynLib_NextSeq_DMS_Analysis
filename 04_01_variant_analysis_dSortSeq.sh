#!/bin/bash

while getopts 'm:M:b:a:"' opt; do
    case "$opt" in

        m)
          min="$OPTARG"
          ;;
        M)
          max="$OPTARG"
          ;;
        b)
          bq="$OPTARG"
          ;;
        a)
          abundance="$OPTARG"
          ;;

    esac
done


for downsample_num in $(seq "$min" "$max")
do

  python3 04_01_variant_analysis_dSortSeq.py -bq "$bq" -a "$abundance" -n "$downsample_num"


done




