#!/bin/bash
#Author: Andrew Nelson andrewnelson@email.arizona.edu
#A bash script to add base pairs to the start and stop coordinates of a genome annotation file (GFF/GFF3/GTF) 
usage(){
    echo ""
    echo "Usage : sh $0 -g gff_file -d delimiter -v base_pairs -o output_file"
    echo ""

cat <<'EOF'

    -g <gff_file>

    -d <delimiter gene/CDS/Gene/mRNA>

    -v <basepairs to modify>

    -o <output_file>

EOF
   exit 0
}

while getopts ":hg:v:d:o:" opt; do
  case $opt in
    g)
      gfftomodify=$OPTARG #File to be modified
    ;;
    d)
      delimiter=$OPTARG #The label provided by the GTF/GFF file denoting the gene or CDS in which you are interested. For instance, most annotation files use gene, but some use Gene or CDS or mRNA.
    ;;
    v)
      basepairtomodify=$OPTARG #This is the number of base pairs that you want to add to the 5' and 3' of your gene. Both negative and positive values can be used here.
    ;;
    o)
      output=$OPTARG # This is output file
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
grep -v "#" $gfftomodify |grep -i "$delimiter	" |awk -v num="$basepairtomodify" '{OFS="\t"} {$4-=num}1' |awk -v num="$basepairtomodify" '{OFS="\t"} {$5+=num}1' > $output
