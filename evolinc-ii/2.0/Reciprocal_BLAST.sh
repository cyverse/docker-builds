#!/bin/bash
#author: Andrew Nelson; andrew.d.l.nelson@gmail.com
#This script is designed to perform reciprocal BLAST of identified lincRNA orthologs in conjunction with Evolinc II
# Usage: 
# sh Reciprocal_BLAST.sh 

usage() {
      echo ""
      echo "Usage : sh $0 -g genome -s ortho_sequences -f ortho_gff -a query_gff -b query_sp -c query_genome"
      echo ""

cat <<'EOF'
     -g    </path/to/genome in FASTA format>

     -s    </path/to/putative ortholog sequences in FASTA format>

     -f    </path/to/gff file of ortholog sequences>

     -a    <query gff in four letter format. ie., Atha for A. thaliana>

     -b    <query species in four letter format>

     -c    <query genome in fasta format>"

     -h    Show this usage information
EOF
    exit 0
}


while getopts ":g:s:f:a:hb:c:" opt; do
  case $opt in
    g)
     Genome=$OPTARG
     ;;
    s)
     Put_ortholog=$OPTARG
     ;;
    f)
     Put_ortholog_gff=$OPTARG
     ;;
    a)
     query_gff=$OPTARG
     ;;
    b)
     query_species=$OPTARG
     ;;
    c)
     query_genome=$OPTARG
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

echo "Starting on $Put_ortholog"
mkdir -p Reciprocal_BLAST_DB
mkdir -p Reciprocal_BLAST_Return

# Formatting the genome
makeblastdb -logfile stderr.out -in $query_genome -dbtype nucl -out Reciprocal_BLAST_DB/$query_genome.blast.out

# Blasting the putative ortholog sequences against the query genome genome to find out the location on the genome.
blastn -logfile stderr.out -query $Put_ortholog -db Reciprocal_BLAST_DB/$query_genome.blast.out -num_threads 4 -penalty -2 -reward 1 -gapopen 5 -gapextend 2 -dust no -word_size 8 -evalue 1e-20 -outfmt "6 qseqid sseqid pident length qlen qstart qend sstart send evalue bitscore" -out Reciprocal_BLAST_Return/$Put_ortholog.Reciprocal.out

# Remove spaces in the blastout files
sed 's/ //g' Reciprocal_BLAST_Return/$Put_ortholog.Reciprocal.out > Reciprocal_BLAST_Return/$Put_ortholog.Reciprocal.stripped.out

# Convert blast result to gff
perl /blast2gff_recip.pl -i Reciprocal_BLAST_Return/$Put_ortholog.Reciprocal.stripped.out
# Change the file format since the gffread requires the chromosome to be first column
awk '{print $2 "\t" $1 "\t" $3 "\t" $4 "\t" $5 "\t" $6 "\t" $7 "\t" $8 "\t" $9 " " $10  " " $11 " " $12}' Reciprocal_BLAST_Return/$Put_ortholog.Reciprocal.stripped.out.gff > temp && mv temp Reciprocal_BLAST_Return/$Put_ortholog.Reciprocal.stripped.out.gff

grep "TBH" Reciprocal_BLAST_Return/$Put_ortholog.Reciprocal.stripped.out.gff >Reciprocal_BLAST_Return/$Put_ortholog.Reciprocal.TBH.out.gff
bedtools sort -i Reciprocal_BLAST_Return/$Put_ortholog.Reciprocal.TBH.out.gff >Reciprocal_BLAST_Return/$Put_ortholog.sorted.gff
cut -f 2 Reciprocal_BLAST_Return/$Put_ortholog.sorted.gff >Reciprocal_BLAST_Return/$Put_ortholog.reciprocal.Unique.IDs.almost.txt
grep "TBH" $query_gff >Reciprocal_BLAST_Return/$query_gff.TBH.gff
bedtools sort -i Reciprocal_BLAST_Return/$query_gff.TBH.gff >Reciprocal_BLAST_Return/$query_gff.TBH.sorted.gff
bedtools closest -a Reciprocal_BLAST_Return/$query_gff.TBH.sorted.gff -b Reciprocal_BLAST_Return/$Put_ortholog.sorted.gff >Reciprocal_BLAST_Return/$Put_ortholog.reciprocal.comparison.gff
cut -f 1,2,4,5,10,11,13,14 Reciprocal_BLAST_Return/$Put_ortholog.reciprocal.comparison.gff >Reciprocal_BLAST_Return/$Put_ortholog.reciprocal.comparison.filtered.gff
#The below step is causing issues and I cannot remember why I inserted it here.
#grep -f Reciprocal_BLAST_Return/$Put_ortholog.reciprocal.Unique.IDs.list.txt Reciprocal_BLAST_Return/$Put_ortholog.reciprocal.comparison.filtered.gff >Reciprocal_BLAST_Return/$Put_ortholog.reciprocal.comparison.TBH.filtered.gff

sed -i 's~.....TCONS~TCONS~g' Reciprocal_BLAST_Return/$Put_ortholog.reciprocal.comparison.filtered.gff
awk -F'\t' '{gsub(/_Known.*/,"",$6); print}' Reciprocal_BLAST_Return/$Put_ortholog.reciprocal.comparison.filtered.gff > temp && mv temp Reciprocal_BLAST_Return/$Put_ortholog.reciprocal.comparison.filtered.gff

###
awk -F'\t' 'substr($2,6) == substr($6,6) {print $2}' Reciprocal_BLAST_Return/$Put_ortholog.reciprocal.comparison.filtered.gff >Reciprocal_BLAST_Return/$Put_ortholog.reciprocal.IDs.list.txt
#awk -F'\t' '{ if (substr($2,6) == substr($6,6)) { print $2} else { pass;} }' Reciprocal_BLAST_Return/$Put_ortholog.reciprocal.comparison.filtered.gff >temp && mv temp Reciprocal_BLAST_Return/$Put_ortholog.reciprocal.IDs.list.txt
uniq Reciprocal_BLAST_Return/$Put_ortholog.reciprocal.IDs.list.txt >Reciprocal_BLAST_Return/$Put_ortholog.reciprocal.Unique.IDs.list.txt
#The line below is being used to generate a list of families (and all their paralogous sequences) to be pulled from the overall FASTA file of sequence homologs in each species.
sed 's~_TBH_1~~g' Reciprocal_BLAST_Return/$Put_ortholog.reciprocal.Unique.IDs.list.txt >Reciprocal_BLAST_Return/$Put_ortholog.reciprocal.final.list.for.pooling.families.txt
#assign subject_species argument
subject_species=${Put_ortholog:0:4}
#sed -i 's~=.*~~g' Reciprocal_BLAST_Return/$Put_ortholog.reciprocal.final.list.for.pooling.families.txt
sed -i 's~^~>'$subject_species'_~g' Reciprocal_BLAST_Return/$Put_ortholog.reciprocal.final.list.for.pooling.families.txt

sleep 5
sed 's~^.....~~g' Reciprocal_BLAST_Return/$Put_ortholog.reciprocal.Unique.IDs.list.txt > Reciprocal_BLAST_Return/$Put_ortholog.reciprocal.Unique.IDs.stripped.list.txt
grep -Ef Reciprocal_BLAST_Return/$Put_ortholog.reciprocal.Unique.IDs.stripped.list.txt Reciprocal_BLAST_Return/$Put_ortholog.reciprocal.Unique.IDs.almost.txt >Reciprocal_BLAST_Return/$Put_ortholog.reciprocal.comparison.TBH.filtered.gff
sort -u Reciprocal_BLAST_Return/$Put_ortholog.reciprocal.comparison.TBH.filtered.gff >Reciprocal_BLAST_Return/$Put_ortholog.reciprocal.Unique.IDs.for.renaming.final.txt

find ./Reciprocal_BLAST_Return -type f ! -name '*final*' -delete
### Error here now!
python /assign_annotation_ortholog.py ../Orthologs/$Put_ortholog Reciprocal_BLAST_Return/$Put_ortholog.reciprocal.Unique.IDs.for.renaming.final.txt ../Orthologs/$Put_ortholog.orthologs.identified.fasta
rm -f ../Orthologs/$Put_ortholog
perl /singleline.pl ../Orthologs/$Put_ortholog.orthologs.identified.fasta >../Orthologs/$Put_ortholog.orthologs.identified.singleline.fasta
grep -A 1 -f Reciprocal_BLAST_Return/$Put_ortholog.reciprocal.final.list.for.pooling.families.txt ../Orthologs/$Put_ortholog.orthologs.identified.singleline.fasta >../Orthologs/$Put_ortholog.orthologs.only.fasta
rm -f ../Orthologs/$Put_ortholog.orthologs.identified.fasta
rm -f ../Orthologs/$Put_ortholog.orthologs.identified.singleline.fasta
find ./Reciprocal_BLAST_Return -type f ! -name '*renaming*' -delete
echo "Finished with $Put_ortholog"
