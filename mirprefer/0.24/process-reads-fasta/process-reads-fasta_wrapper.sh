#!/bin/bash

usage() {
      echo ""
      echo "Usage : sh $0 -s SAMPLE_LIST -f FASTA_FILE1 [-f FASTA_FILE2 ....] -o OUTPUT_FOLDER"
      echo ""

cat <<'EOF'

  -s </path/to/sample_list file> [MANDATORY]
  -f </path/to/fasta file (s)> [Mandatory]
  -o </path/to/output folder> [Mandatory]

EOF
    exit 0
}


while getopts ":s:hf:o:" opt; do
  case $opt in
    f)
     fastq_file+=("$OPTARG")
      ;;
    h)
     usage
     exit 1
      ;;    
    s)
     sample_list=$OPTARG
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

# Convert fastq to fasta
for f in "${fastq_file[@]}"; do
    new=$(basename $f ".fastq")
    cat $f | paste - - - - | sed 's/^@/>/g' | cut -f1-2 | tr '\t' '\n' | sed 's/ .*//g' > ${new}.fasta
done

# Process reads fasta
python /miR-PREFeR-0.24/scripts/process-reads-fasta.py "$sample_list" *.fasta

# Remove fasta files and also move output to a folder
rm *.fasta
mv *processed $output

