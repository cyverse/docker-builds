#!/bin/bash
# A script necessary to merge GTF entries that overlap.
# Author: Andrew Nelson
# Date: 10/23/17

usage() {
      echo ""
      echo "Usage : sh $0 -c cuffcompare output file -o filtered output file" 
      echo ""

cat <<'EOF'
  -c </path/to/cuffcompare output file>
  -o <output file name>
EOF
    exit 0
}  
while getopts ":c:o:" opt; do
  case $opt in
    c)
     comparefile=$OPTARG
     ;;
    o)
     outputfile=$OPTARG
     ;;  
esac
done
  
grep -v -w '\W+\W' $comparefile | grep -v -w '\W-\W' >unstranded.gtf
grep -w '\W+\W' $comparefile >positive_stranded.gtf
grep -w '\W-\W' $comparefile >negative_stranded.gtf

intersectBed -a unstranded.gtf -b positive_stranded.gtf -u >positive_unstranded_overlap.gtf
cut -f 9 positive_unstranded_overlap.gtf | awk -F " " '{print $2}' | sort | uniq | sed 's~;~~g' >unstranded_to_remove.txt
grep -vFf unstranded_to_remove.txt unstranded.gtf >unstranded_clean.gtf
awk -F'\t' '{$7="+"} 1' OFS='\t' unstranded_clean.gtf >unstranded_now_stranded.gtf
cat unstranded_clean.gtf positive_stranded.gtf >positive_plus_unstranded.gtf
sortBed -i positive_plus_unstranded.gtf >positive_plus_unstranded_sorted.gtf
intersectBed -a negative_stranded.gtf -b positive_plus_unstranded_sorted.gtf >negative_positive_overlap.gtf
intersectBed -a positive_stranded.gtf -b negative_stranded.gtf -u >positive_stranded_for_renaming.gtf
cut -f 9 negative_positive_overlap.gtf | awk -F " " '{print $2}' | sort | uniq | sed 's~;~~g' >negative_to_remove.txt
cut -f 9 positive_stranded_for_renaming.gtf | awk -F " " '{print $2}' | sort | uniq | sed 's~;~~g'| tr -d '"' >positive_stranded_list.txt
sed 's~XLOC_~XLOC.uncertain.strand.support_~g' positive_stranded_list.txt >positive_stranded_list_second_column.txt
paste positive_stranded_list.txt positive_stranded_list_second_column.txt >renaming_list.txt
grep -vFf negative_to_remove.txt negative_stranded.gtf >negative_stranded_clean.gtf
cat negative_stranded_clean.gtf positive_plus_unstranded_sorted.gtf >cuffmerge_collapsed.gtf
sortBed -i cuffmerge_collapsed.gtf >$outputfile
# rm positive_* negative_* unstranded* cuffmerge_collapsed.gtf