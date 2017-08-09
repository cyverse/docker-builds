#!/bin/bash

usage() {
      echo ""
      echo "Usage : sh $0 -i input_file"
      echo ""

cat <<'EOF'
  -i </path/to/bam_file>

  -h Show this usage information

EOF
    exit 0
}

while getopts ":hi:" opt; do
  case $opt in
    i)
     input=$OPTARG
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

for file in $input/*; 
   do 
   var="$(samtools view -F 4 -c $file)"; 
   var1="$(samtools view -f 4 -c $file)";
   var3="$(samtools view -c $file)";
   echo -e $file '\t' $var '\t' $var1 '\t' $var3 >> bam_stats.txt; 
done

awk 'BEGIN {print "Sample_name\tMapped_reads\tUnmapped_reads\tTotal_reads"} {print}' bam_stats.txt > temp && mv temp bam_stats.txt
