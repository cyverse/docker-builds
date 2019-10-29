#!/bin/bash

set -x

usage() {
      echo ""
      echo "Usage : sh $0 -d <Maker-out-log-file>"
      echo ""
}

while getopts ":hd:" opt; do
  case $opt in
    d)
    maker_out=$OPTARG # Reference genome file
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

# gff3_merge
gff3_merge -d $maker_out/*index.log -o my_assembly.all.gff

# maker2zff conversion
maker2zff my_assembly.all.gff

# fathom
fathom genome.ann genome.dna -categorize 1000
fathom uni.ann uni.dna -export 1000 -plus 

# forge
forge export.ann export.dna

# hmm-assembler
hmm-assembler.pl my_genome . > my_genome.hmm
