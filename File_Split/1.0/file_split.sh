#!/bin/bash

usage() {
      echo ""
      echo "Usage : sh $0 -i <input file> -f <target folder> -l <number of lines> -p <prefix>"
      echo ""
	  cat <<'EOF'
	  ###### Command line options ##########
  -i <input file>
  -f <folder location (target folder) in which to place the files>
  -l <number of lines per output file>
  -p <prefix of the file>
EOF
    exit 0
}
while getopts ":hf:i:l:p:" opt; do
  case $opt in
    i)
    inputfile=$OPTARG #Input directory where files are currently located
     ;;
    f)
    foldername=$OPTARG # Target directory
     ;;
    l)
    lines=$OPTARG # Number of lines per output file
     ;;
    p)
    prefix=$OPTARG # Prefix of the files e.g: test_ (Default is x)
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

# Run split command
split $inputfile -l $lines $prefix

# Move the files into target directory
find . -maxdepth 1 \( ! -type d \) -exec mv {} ./$foldername \;
