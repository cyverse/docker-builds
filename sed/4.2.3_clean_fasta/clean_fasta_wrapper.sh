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
     input+=("$OPTARG")
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

mkdir $output

for f in "${input[@]}"; do
    extension=${f#*.}
    if [ "$extension" == "fasta" ]; then
        filename=$(basename "$f" ".fasta")
        sed 's/.*|/>/g' $f > out.txt && mv out.txt $output/"$filename".cleaned."$extension"
    elif [ "$extension" == "fa" ]; then
        filename=$(basename "$f" ".fa")
        sed 's/.*|/>/g' $f > out.txt && mv out.txt $output/"$filename".cleaned."$extension"
   else
	echo The extension - "$extension" for the file "$f" needs to be either fasta or fa
   fi

done
