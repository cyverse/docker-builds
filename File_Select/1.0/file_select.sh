#!/bin/bash

usage() {
      echo ""
      echo "Usage : sh $0 -i <input directory> -f <target folder> -a <pattern1> -a <pattern2>"
      echo ""
	  cat <<'EOF'
	  ###### Command line options ##########
  -i <input directory>
  -f <folder location (target folder) in which to place the files>
  -a <one or more patterns eg. *combined.gtf or *sorted.bam>
EOF
    exit 0
}
while getopts ":hf:i:a:" opt; do
  case $opt in
    i)
    input_directory=$OPTARG #Input directory where files are currently located
     ;;
    f)
    foldername=$OPTARG # Target directory
     ;;
    a)
    filetype+=("$OPTARG") # first pattern to search for, e.g., *combined.gtf, be sure to include the asterisk
     ;;
    h)
    usage
     exit 1
      ;;
    \?)
      echo "Invalid option: -$OPTARG" >&2
      exit 1
      ;;
    :)
      echo "Option -$OPTARG requires an argument." >&2
      exit 1
      ;;
  esac
done

# Create the targe directory
mkdir -p $foldername

# Loop over the patterns in the input directory and move the files into target directory
for type in "${filetype[@]}"
do
        find $input_directory/. -name "$type" -exec mv {} ./$foldername \;
done
