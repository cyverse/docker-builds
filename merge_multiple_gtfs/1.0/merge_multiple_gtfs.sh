#!/bin/bash
#Script to merge multiple Evolinc updated GTF files into one

#find all Evolinc output GTF files
# As input this script requires the folder where all of the updated GTFs are located (they can be in subfolders)
## in addition, this script requires the GTF that was used as input for RMTA - not Evolinc.
### This script also requires the original genome (.fasta) used in Evolinc or RMTA
#### Prerequisite programs this script needs are bedtools and cufflinks (cuffmerge)

#Assign input RMTA-GTF to $originalgtf
#Assign input genome to $genome

usage() {
      echo ""
      echo "Usage : sh $0 -g reference genome -o reference gtf -d directory"
      echo ""

cat <<'EOF'
  -g </path/to/genome file>

  -o </path/to/refefence gtf>

  -d </parent/directory/of/gtfs>

EOF
    exit 0
}  
while getopts ":g:o:d:" opt; do
  case $opt in
    g)
     genome=$OPTARG
      ;;  
    o)
     gtf=$OPTARG
      ;;
    d)
     directory=$OPTARG
      ;;
    h)
    usage
    exit 1
      ;;  
esac
done

# Find the updated lincRNA gtf files

find $directory -name "*lincRNA.updated.gtf" -print > list_of_gtfs.txt
cat list_of_gtfs.txt | while read line 
do 
  grep 'lincRNA' $line > temporary.file
  lincRNA=$(echo $line | sed 's~\.~~g' | sed 's~\/~~g')
  echo "assigned new variable"
  mv temporary.file $lincRNA.lincRNA.for.merging
  grep 'uncertain' $line > $lincRNA.unstranded.only
done

# lincRNA

find . -name "*lincRNA.for.merging*" -print > list_of_lincRNAs.txt

cuffmerge -s $genome list_of_lincRNAs.txt

find -name "merged.gtf" -exec mv {} ./all_lincRNAs_merged.gtf \;

sortBed -i all_lincRNAs_merged.gtf > all_lincRNAs_merged_sorted.gtf

# unstranded

find . -name "*unstranded*" -print > list_of_unstranded.txt

cuffmerge -s $genome list_of_unstranded.txt

find -name "merged.gtf" -exec mv {} ./unstranded_lincRNAs_merged.gtf \;

sortBed -i unstranded_lincRNAs_merged.gtf > unstranded_lincRNAs_merged_sorted.gtf

# Renaming

intersectBed -a all_lincRNAs_merged_sorted.gtf -b unstranded_lincRNAs_merged_sorted.gtf > lincRNAs_for_renaming.gtf

cut -f 9 lincRNAs_for_renaming.gtf | awk -F " " '{print $2}' | sort | uniq | sed 's~;~~g'| tr -d '"' > lincRNA_first_column.txt

sed  's~XLOC_~XLOC_uncertain_strand_support_~g' lincRNA_first_column.txt > lincRNA_second_column.txt

paste lincRNA_first_column.txt lincRNA_second_column.txt > lincRNA_renaming_list.txt

while read -r pattern replacement; do
   sed -i "s/$pattern/$replacement/" all_lincRNAs_merged_sorted.gtf
done < lincRNA_renaming_list.txt

sed 's~Cufflinks~Evolinc~g' all_lincRNAs_merged_sorted.gtf > Updated_lincRNAs.gtf

# Extract the fasta sequences of the lincRNA

gffread -w Updated_lincRNAs.fa -g $genome Updated_lincRNAs.gtf

#Clean up lincRNA IDs

sed -i 's~ gene=~~g' Updated_lincRNAs.fa
sed -i 's~TCONS_~lincRNA.ID.~g' Updated_lincRNAs.fa
sed -i 's~XLOC_~.locus.ID.~g' Updated_lincRNAs.fa

# Update the reference gtf file

cat Updated_lincRNAs.gtf $gtf > Updated.gtf

sortBed -i Updated.gtf > Final_updated.gtf

# Clean up

rm list_of_lincRNAs.txt list_of_gtfs.txt *for.merging all_lincRNAs_merged.gtf all_lincRNAs_merged_sorted.gtf list_of_unstranded.txt unstranded_lincRNAs_merged.gtf unstranded_lincRNAs_merged_sorted.gtf lincRNAs_for_renaming.gtf lincRNA_first_column.txt 
rm *unstranded.only lincRNA_second_column.txt lincRNA_renaming_list.txt Updated_lincRNAs.gtf Updated.gtf -r merged_asm
