#!/bin/bash

usage() {
       echo ""
      echo "Usage : sh $0 -f FASTA_FILE1 [-f FASTA_FILE2 ....] -r REFERENCE_GENOME -t NUMBER_OF_THREADS [-c] -o OUTPUT_FOLDER"
      echo ""

cat <<'EOF'

  -f </path/to/fasta file (s)> [Mandatory]
  -a </path/to/reference genome> [Mandtory]
  -t <number of threads> [Mandatory]
  -c <Filter SAM file to remove unmapped reads or not> [Optional]
  -o </path/to/output folder> [Mandatory]

EOF
    exit 0
}

filt=0

while getopts ":r:hf:ct:o:" opt; do
  case $opt in
    f)
     fasta_file+=("$OPTARG")
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
    c)
     filt=$OPTARG
      ;;
    o)
     output=$OPTARG
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

# Make directory
mkdir $output

# Indexing the reference genome
bowtie2-build $ref_genome reference_index

for f in "${fasta_file[@]}"; do
    # Mapping the reads to the reference genome
    bowtie2 -p $threads -f -x reference_index -U $f -S ${f}.sam
    # Filtering the sam file to remove ummapped reads
    if [ "$filt" != 0 ]; then
       samtools view -Sh -F 4 -o ${f}.filt.sam ${f}.sam
       mv ${f}.sam ${f}.filt.sam $output
    else 
       mv ${f}.sam $output
    fi
done

mv reference_index* $output
