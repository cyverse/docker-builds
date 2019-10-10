#!/bin/bash

usage() {
       echo ""
      echo "Usage : sh $0 -f FASTA_FILE -r REFERENCE_GENOME -t NUMBER_OF_THREADS"
      echo ""

cat <<'EOF'

  -f </path/to/fasta file>
  -a </path/to/reference genome>
  -t <number of threads>

EOF
    exit 0
}

while getopts ":r:hf:t:" opt; do
  case $opt in
    f)
     fasta_file=$OPTARG
      ;;
    h)
     usage
     exit 1
      ;;    
    r)
     ref_genome=$OPTARG
      ;;
    t)
     threads=$OPTARG
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

bowtie2-build $ref_genome reference_index
bowtie2 -p $threads -f -x reference_index -U $fasta_file -S ${fasta_file}.sam &
samtools view -Sh -F 4 -o ${fasta_file}.filt.sam ${fasta_file}.sam
