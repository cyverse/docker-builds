#!/bin/bash

usage() {
      echo ""
      echo "Usage : sh $0 -i <input SRA_ID file> -o <output folder> -p <number of threads>"
      echo ""

cat <<'EOF'
  -i </path/to/input sraid file>
  -o </path/to/output folder>
  -p <number of threads>
EOF
    exit 0
}

while getopts ":i:o:p:hn:x:" opt; do
  case $opt in
    i)
     input=$OPTARG
      ;;
    o)
     output=$OPTARG
      ;;
    p)
     threads=$OPTARG
      ;;
    # n)
    #  minspotid=$OPTARG
    #   ;;
    # x)
    #  maxspotid=$OPTARG
    #   ;;
    h)
     usage
     exit 1
      ;;    
  esac
done

mkdir $output

while read f; do

  mkdir $f.output
  parallel-fastq-dump -s $f -t $threads -O $f.output --split-files --gzip
  mv $f.output/* $output
  rm -r $f.output
done < "$input"



