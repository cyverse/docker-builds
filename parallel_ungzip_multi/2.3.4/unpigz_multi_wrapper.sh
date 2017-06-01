#!/bin/bash

usage() {
      echo ""
      echo "Usage : sh $0 [-f <faster>] [-b <best>] [-p <number of threads>] -i <input file> -o <output folder>"
      echo ""

cat <<'EOF'

  -i </path/to/input file>
 
  -f <compress faster>

  -b <compress better>

  -p <number of threads>
 
  -o </path/to/output file>

EOF
    exit 0
}

fast=0
best=0

while getopts ":i:f:b:p:ho:" opt; do
  case $opt in
    i)
     input+=("$OPTARG")
      ;;
    f)
     fast=$OPTARG
      ;;
    b)
     best=$OPTARG
      ;;
    p)
     proc=$OPTARG
      ;;
    o)
     output=$OPTARG
      ;;
    h)
     usage
     exit 1
      ;;    
  esac
done

mkdir $output

if [ ! -z $input ] && [ ! -z $output ] && [ $fast == 0 ] && [ $best == 0 ]; then
	for f in "${input[@]}"; do
	     file=$(basename "$f" ".gz")
	     unpigz -p $proc $f
	     mv $file $output 
	done
elif [ ! -z $input ] && [ ! -z $output ] && [ $fast != 0 ] && [ $best == 0 ]; then
        for f in "${input[@]}"; do
	     file=$(basename "$f" ".gz")
	     unpigz --fast -p $proc $f
	     mv $file $output 
        done
elif [ ! -z $input ] && [ ! -z $output ] && [ $fast == 0 ] && [ $best != 0 ]; then
        for f in "${input[@]}"; do
	     file=$(basename "$f" ".gz")
             unpigz --best -p $proc $f
	     mv $file $output 
        done
fi
