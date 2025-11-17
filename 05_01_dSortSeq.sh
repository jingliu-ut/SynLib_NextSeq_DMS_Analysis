#!/bin/bash

while getopts 'c:n:' opt; do
    case "$opt" in
        c)
          condition="$OPTARG"
          ;;
        n)
          downsample_num="$OPTARG"
    esac
done


  if [ "$condition" != "all" ]; then

    input_dir="/media/scratch/jingscratch/SynLib_NextSeq_DMS/04_01_dSortSeq_inputs/downsample_${downsample_num}"
    output_dir="/media/scratch/jingscratch/SynLib_NextSeq_DMS/05_01_dSortSeq_outputs/downsample_${downsample_num}/${condition}"
    mkdir -p "$output_dir"

    mixing_coefficient="${input_dir}/${condition}_mixing_coefficient.csv"
    binned_distribution="${input_dir}/${condition}_binned_distribution.csv"
    boundaries="${input_dir}/${condition}_boundaries.csv"
    overall_distribution="${input_dir}/${condition}_overall_distribution.csv"

    python 05_dSortSeq.py -p "$mixing_coefficient" -f "$binned_distribution" -d "$overall_distribution" -b "$boundaries" -o "$output_dir" -c "$condition"

  fi

  # Check if -c is "all"
  if [ "$condition" == "all" ]; then

      conditions=("dark" "light" "ymScarlet")

      for condition in "${conditions[@]}"; do

          input_dir="/media/scratch/jingscratch/SynLib_NextSeq_DMS/04_01_dSortSeq_inputs/downsample_${downsample_num}"
          output_dir="/media/scratch/jingscratch/SynLib_NextSeq_DMS/05_01_dSortSeq_outputs/downsample_${downsample_num}/${condition}"
          mkdir -p "$output_dir"

          mixing_coefficient="${input_dir}/${condition}_mixing_coefficient.csv"
          binned_distribution="${input_dir}/${condition}_binned_distribution.csv"
          boundaries="${input_dir}/${condition}_boundaries.csv"
          overall_distribution="${input_dir}/${condition}_overall_distribution.csv"

          python 05_dSortSeq.py -p "$mixing_coefficient" -f "$binned_distribution" -d "$overall_distribution" -b "$boundaries" -o "$output_dir" -c "$condition"

          echo "Waiting to finish..."
          wait

          for i in $(seq 1 30); do
              echo "Sleep $i sec..."
              sleep 1
          done
      done

      condition="all"
  fi
