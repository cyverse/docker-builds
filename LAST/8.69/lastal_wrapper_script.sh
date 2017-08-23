#!/bin/bash

set -x
set -u

while getopts ":hP:m:E:i:q:" opt; do
  case $opt in
  	P)
	threads=$OPTARG
	;;
	m)
	sensitivity=$OPTARG
	;;
	E)
	evalue=$OPTARG
	;;
	i)
	index=$OPTARG
	;;
	q)
	query+=("$OPTARG")
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
   
# Index
fil=$(ls $index/* | head -n 1)
filename=$(basename "$fil")
filename="${filename%.*}"
cp $index/* .

# Input file
filename2=$(basename "$query")
filename2="${filename2%.*}"

# LASTAL 
lastal -P $threads -E $evalue -m $sensitivity $filename $query | last-split -m1 > $filename2.out.maf
maf-swap $filename2.out.maf | last-split -m1 > $filename2.out2.maf

# Clean up
rm $filename.*

