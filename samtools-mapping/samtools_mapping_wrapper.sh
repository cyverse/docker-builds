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
	echo $file >> output1; 
	samtools view -F 4 $file | wc -l >> output1;
	tr '\n' '\t' < output1 | sed -r 's/([0-9]+\t)/\1\n/' > output2
	done
rm output1
