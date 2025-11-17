downsample_root="/media/scratch/jingscratch/SynLib_NextSeq_DMS/02_01_downsample_outputs"

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

# Check if downsample_none exists in 03_01_vcf_files
if [ ! -d "/media/scratch/jingscratch/SynLib_NextSeq_DMS/03_01_vcf_files/downsample_none" ]; then
  echo "downsample_none does not exist in 03_01_vcf_files."

  echo "Running 03_variant_calling.py for downsample_none"

  input_dir="$downsample_root/downsample_none"
  mkdir -p "/media/scratch/jingscratch/SynLib_NextSeq_DMS/03_01_vcf_files/downsample_none"


  process_files() {
  local start=$1
  local end=$2

  for file in "$input_dir"/*.fastq.gz
  do
    filename=${file##*/}
    sample=${filename%%_*}
    if (( sample >= start && sample <= end )); then
      echo "Running 03_variant_calling.py for $file"

      python3 03_variant_calling.py -n -1 -i "$file"
    fi
  done
  }

  process_files 1 4 &
  process_files 8 13 &
  process_files 16 20 &
  process_files 25 25 &

  wait

fi



for i in $(seq "$min" "$max")
do
  # Run 03_variant_calling.py
  echo "Running 03_variant_calling.py"
  input_dir="$downsample_root/downsample_$i"
  output_dir="/media/scratch/jingscratch/SynLib_NextSeq_DMS/03_01_vcf_files/downsample_$i"
  mkdir -p "$output_dir"

  process_files() {
  local start=$1
  local end=$2

  for file in "$input_dir"/*.fastq.gz
  do
    filename=${file##*/}
    sample=${filename%%_*}
    if (( sample >= start && sample <= end )); then
      echo "Running 03_variant_calling.py for $file"

      python3 03_variant_calling.py -n "$i" -i "$file"
    fi
  done
  }

  process_files 2 5 &
  process_files 6 9 &
  process_files 10 12 &
  process_files 14 17 &
  process_files 18 21 &
  process_files 22 24 &

  wait
done
wait