#!/bin/bash

usage() {
      echo ""
      echo "Usage : sh $0 -i <input file> -o <output file>"
      echo ""

cat <<'EOF'

  -i </path/to/input filee>

  -o </path/to/output file>

EOF
    exit 0
}

while getopts ":i:ho:" opt; do
  case $opt in
    i)
     input=$OPTARG
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

sed 's/.*|/>/g' $input > $output
