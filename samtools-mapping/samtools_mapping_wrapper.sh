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

#!/bin/bash
for file in test_data/*.bam; 
do 
var="$(samtools view -F 4 $file | wc -l)"; 
var1="$(samtools view -f 4 $file | wc -l)";
var3="$(samtools view $file | wc -l)";
echo -e $file '\t' $var '\t' $var1 '\t' $var3 >> output1.txt; 
done

awk 'BEGIN {print "Sample_name\tMapped_reads\tUnmapped_reads\tTotal_reads"} {print}' output1.txt > output.txt
rm output1.txt
