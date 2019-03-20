#!/bin/bash

set -x
set -e

while getopts ":hi:t:" opt; do
  case $opt in
    i)
    input=$OPTARG # Input fasta file
     ;;
    t)
    threads=$OPTARG # Number of cpus
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

# Total residues
total=$(esl-seqstat -a $input | grep "Total" | tr ':' '\n' | grep [0-9])

# Converting bp into Mbp
final=$(echo "$total * 2 / 1000000" | bc -l)

# Extract
cmpress /Rfam.cm

# Extracting the prefix
new=$(basename $input | cut -d . -f 1)

# Running cmscan
cmscan -Z $final --cut_ga --rfam --nohmmonly --tblout "$new".tblout --fmt 2 --cpu $threads --clanin /Rfam.clanin /Rfam.cm $input > "$new".cmscan

