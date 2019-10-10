#!/bin/bash

usage() {
      echo ""
      echo "Usage : sh $0 -i <input file> -d <de file> -o <output file>"
      echo ""

cat <<'EOF'

  -i </path/to/input custom file>
 
  -d </path/to/de genome file>

  -o </path/to/output file>

EOF
    exit 0
}

while getopts ":i:d:ho:" opt; do
  case $opt in
    i)
     input+=("$OPTARG")
      ;;
    d)
     de+=($OPTARG)
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

if [ ! -z $input ]; then
	for f in "${input[@]}"; do
    	    extension=${f#*.}
    	    if [ "$extension" == "fasta" ]; then
        	filename=$(basename "$f" ".fasta")
        	sed 's/.*|/>/g' $f > out.txt && mv out.txt $output/"$filename".cleaned."$extension"
            elif [ "$extension" == "fa" ]; then
        	filename=$(basename "$f" ".fa")
           	sed 's/.*|/>/g' $f > out.txt && mv out.txt $output/"$filename".cleaned."$extension"
    	    elif [ "$extension" == "fas" ]; then
        	filename=$(basename "$f" ".fas")
        	sed 's/.*|/>/g' $f > out.txt && mv out.txt $output/"$filename".cleaned."$extension"
    	    else
		echo The extension - "$extension" for the file "$f" needs to be either fasta or fas or fa 1>&2
		rm -r $output
		exit 64
   	    fi
	done
elif [ ! -z $de ]; then
        cp $de .
        file=$(basename "$de")
	extension=${file#*.}
        if [ "$extension" == "fasta" ]; then
           filename=$(basename "$file" ".fasta")
           sed 's/.*|/>/g' "$file" > out.txt && mv out.txt $output/"$filename".cleaned."$extension"
        elif [ "$extension" == "fa" ]; then
           filename=$(basename "$file" ".fa")
           sed 's/.*|/>/g' "$file" > out.txt && mv out.txt $output/"$filename".cleaned."$extension"
        elif [ "$extension" == "fas" ]; then
           filename=$(basename "$file" ".fas")
           sed 's/.*|/>/g' "$file" > out.txt && mv out.txt $output/"$filename".cleaned."$extension"
        else
           echo The extension - "$extension" for the file "$f" needs to be either fasta or fas or fa 1>&2
           rm -r $output
           exit 64
        fi
	rm "$file"  
fi

