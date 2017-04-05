#!/bin/bash
# Upendra Kumar Devisetty
# Script to process cuffcompare output file to generate long non-coding RNA

usage() {
      echo ""
      echo "Usage : sh $0 -c cuffcompare -g genome -u user_gff -r gff -o output -n threads [-b TE_RNA] [-t CAGE_RNA] [-x Known_lincRNA]"
      echo ""

cat <<'EOF'

  -c </path/to/cuffcompare output file>

  -g </path/to/reference genome file>

  -u </path/to/user reference annotation file>

  -r </path/to/reference annotation file>

  -b </path/to/Transposable Elements file>

  -n <number of threads>

  -o </path/to/output file>

  -t </path/to/CAGE RNA file>
  
  -x </path/to/Known lincRNA file>

  -h Show this usage information

EOF
    exit 0
}

while getopts ":b:c:g:hr:t:x:o:n:u:" opt; do
  case $opt in
    b)
     blastfile=$OPTARG
      ;;
    c)
     comparefile=$OPTARG
      ;;
    h)
     usage
     exit 1
      ;;    
    g)
     referencegenome=$OPTARG
      ;;
    u)
     user_referencegff=$OPTARG
      ;;
    r)
     referencegff=$OPTARG
      ;;
    n)
     threads=$OPTARG
      ;;  
    t)
     cagefile=$OPTARG
      ;;
    x)
     knownlinc=$OPTARG
      ;;
    o)
     output=$OPTARG
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

echo "BEGIN!"
echo `date`

START_TIME=$SECONDS

# Throw an error if the genomes are coming from NCBI that adds "|" in the headers
if ( grep -q ".*|" $referencegenome ); then
   echo "Your genome file header have pipe characters. Please remove/replace them before proceeding" 1>&2
   exit 64
fi

# Create a directory to move all the output files
mkdir $output

# Fix the gff files if they are coming for coge
grep -v "#" $referencegff |grep -iv "chromosome    " | grep -iv "region    " | grep -iv "scaffold    " > temp && mv temp $referencegff
grep -v "#" $user_referencegff |grep -iv "chromosome    " | grep -iv "region    " | grep -iv "scaffold    " > temp && mv temp $user_referencegff

# Fixing the cuffcompare files as well
sed 's~^~>~g' $comparefile | sed 's~^>0*~>~g' | sed 's~^>Chr0*~>~g' | sed 's~^>Scaffold0*~>~g' | sed 's~^>~~g' > comparefile.gtf
sed 's~^>0*~>~g' $referencegenome | sed 's~^>Chr0*~>~g' | sed 's~^>Scaffold0*~>~g' > referencegenome.fa
 
# STEP 1:
START_TIME_1=$SECONDS

# Extracting classcode u transcripts, making fasta file, removing transcripts > 200 and selecting protein coding transcripts and modify the header to generate genes
grep '"u"' comparefile.gtf | gffread -w transcripts.u.fa -g referencegenome.fa - && python /evolinc_docker/get_gene_length_filter.py transcripts.u.fa putative_intergenic.genes.fa && sed 's/ .*//' putative_intergenic.genes.fa | sed -ne 's/>//p' > putative_intergenic.genes
grep '"x"' comparefile.gtf | gffread -w transcripts.x.fa -g referencegenome.fa - && python /evolinc_docker/get_gene_length_filter.py transcripts.x.fa transcripts.x.filter.fa && sed 's/ .*//' transcripts.x.filter.fa | sed -ne 's/>//p' > transcripts.x.filter.fa.genes 
grep '"s"' comparefile.gtf | gffread -w transcripts.s.fa -g referencegenome.fa - && python /evolinc_docker/get_gene_length_filter.py transcripts.s.fa transcripts.s.filter.fa && sed 's/ .*//' transcripts.s.filter.fa | sed -ne 's/>//p' > transcripts.s.filter.fa.genes
grep '"o"' comparefile.gtf | gffread -w transcripts.o.fa -g referencegenome.fa - && python /evolinc_docker/get_gene_length_filter.py transcripts.o.fa transcripts.o.filter.fa && sed 's/ .*//' transcripts.o.filter.fa | sed -ne 's/>//p' > transcripts.o.filter.fa.genes
grep '"e"' comparefile.gtf | gffread -w transcripts.e.fa -g referencegenome.fa - && python /evolinc_docker/get_gene_length_filter.py transcripts.e.fa transcripts.e.filter.fa && sed 's/ .*//' transcripts.e.filter.fa | sed -ne 's/>//p' > transcripts.e.filter.fa.genes
grep '"i"' comparefile.gtf | gffread -w transcripts.i.fa -g referencegenome.fa - && python /evolinc_docker/get_gene_length_filter.py transcripts.i.fa transcripts.i.filter.fa && sed 's/ .*//' transcripts.i.filter.fa | sed -ne 's/>//p' > transcripts.i.filter.fa.genes
rm referencegenome.fa

# Concatenate all the genes id's from filtered files
cat transcripts.*.filter.fa.genes > transcripts.all.overlapping.filter.genes

# Make new directory
mkdir transcripts_u_filter.fa.transdecoder_dir

# Move the transcript files to this directory transcripts_u_filter.fa.transdecoder_dir
mv transcripts.* transcripts_u_filter.fa.transdecoder_dir/
mv putative_intergenic.* transcripts_u_filter.fa.transdecoder_dir/
# Change the directory
cd transcripts_u_filter.fa.transdecoder_dir

# Run transdecoder now
for file in *filter.fa; do TransDecoder.LongOrfs -t $file; done

# This groups all the longest_orfs.cds and all the longest_orf.pep files into one, in the transdecoder file.
find . -type f -name longest_orfs.cds -exec cat '{}' \; | cat > longest_orfs_cat.cds 
find . -type f -name longest_orfs.pep -exec cat '{}' \; | cat > longest_orfs_cat.pep 

# Blasting the transcripts to uniprot db
blastp -query longest_orfs_cat.pep -db /evolinc_docker/uniprot_sprot.fa -max_target_seqs 1 -outfmt 6 -evalue 1e-5 -num_threads $threads > longest_orfs_cat.pep.blastp

# Genes in the protein coding genes
sed 's/|.*//' longest_orfs_cat.cds | sed -ne 's/>//p' | uniq > longest_orfs.cds.genes

cut -f1 longest_orfs_cat.pep.blastp | cut -d '|' -f 1 | uniq > longest_orfs_cat.pep.blastp.genes

cat longest_orfs.cds.genes longest_orfs_cat.pep.blastp.genes | sort -u > longest_orfs_cat.cds.pep.blastp.genes

# Remove these protein coding genes from the filter file
grep -v -F -f longest_orfs_cat.cds.pep.blastp.genes transcripts.all.overlapping.filter.genes > transcripts.all.overlapping.filter.not.genes  #I added the -F here, it speeds things up quite a bit as it is searching for exact strings.
sed 's/^/>/' transcripts.all.overlapping.filter.not.genes > temp && mv temp transcripts.all.overlapping.filter.not.genes # changed name here 

# Extract fasta file
cat transcripts.*.filter.fa > transcripts.all.overlapping.filter.fa
python /evolinc_docker/extract_sequences.py transcripts.all.overlapping.filter.not.genes transcripts.all.overlapping.filter.fa transcripts.all.overlapping.filter.not.genes.fa 
sed 's/ /./' transcripts.all.overlapping.filter.not.genes.fa > temp && mv temp transcripts.all.overlapping.filter.not.genes.fa

#clean up the folder for the next round of transdecoder
rm longest_orf*

###repeat transdecoder on just the novel transcripts

for file in putative_intergenic.genes.fa; do TransDecoder.LongOrfs -t $file; done

# This groups all the longest_orfs.cds and all the longest_orf.pep files into one, in the transdecoder file.
find . -type f -name longest_orfs.cds -exec cat '{}' \; | cat > longest_orfs_cat.cds 
find . -type f -name longest_orfs.pep -exec cat '{}' \; | cat > longest_orfs_cat.pep 

# Blasting the transcripts to uniprot db
blastp -query longest_orfs_cat.pep -db /evolinc_docker/uniprot_sprot.fa -max_target_seqs 1 -outfmt 6 -evalue 1e-5 -num_threads $threads > longest_orfs_cat.pep.blastp

# Genes in the protein coding genes
sed 's/|.*//' longest_orfs_cat.cds | sed -ne 's/>//p' | uniq > longest_orfs.cds.genes

cut -f1 longest_orfs_cat.pep.blastp | cut -d '|' -f 1 | uniq > longest_orfs_cat.pep.blastp.genes

cat longest_orfs.cds.genes longest_orfs_cat.pep.blastp.genes | sort -u > longest_orfs_cat.cds.pep.blastp.genes

# Remove these protein coding genes from the filter file
grep -v -F -f longest_orfs_cat.cds.pep.blastp.genes putative_intergenic.genes > putative_intergenic.genes.not.genes  #I added the -F here, it speeds things up quite a bit as it is searching for exact strings.
sed 's/^/>/' putative_intergenic.genes.not.genes > temp && mv temp putative_intergenic.genes.not.genes # changed name here 

# Extract fasta file
python /evolinc_docker/extract_sequences.py putative_intergenic.genes.not.genes putative_intergenic.genes.fa putative_intergenic.genes.not.genes.fa 
sed 's/ /./' putative_intergenic.genes.not.genes.fa > temp && mv temp putative_intergenic.genes.not.genes.fa


# Blast the putative intergenic lincRNA fasta file to TE RNA db
if [ ! -z $blastfile ]; then
     makeblastdb -in ../$blastfile -dbtype nucl -out ../$blastfile.blast.out &&
     blastn -query putative_intergenic.genes.not.genes.fa -db ../$blastfile.blast.out -out putative_intergenic.genes.not.genes.fa.blast.out -outfmt 6 -num_threads $threads # no blast hits here
else
    touch putative_intergenic.genes.not.genes.fa.blast.out
fi

# Filter the output to select the best transcript based on e-value and bit-score
python /evolinc_docker/filter_sequences.py putative_intergenic.genes.not.genes.fa.blast.out putative_intergenic.genes.not.genes.fa.blast.out.filtered

# Modify the header in the fasta file to extract header only
grep ">" putative_intergenic.genes.not.genes.fa | sed 's/>//' > putative_intergenic.genes.not.genes_only

# Now remove the blast hits from the fasta file
python /evolinc_docker/fasta_remove.py putative_intergenic.genes.not.genes.fa.blast.out.filtered putative_intergenic.genes.not.genes_only lincRNA.genes

# Modify the fasta header to include ">", generating a new file so as to keep the lincRNA.genes file intact for later use.
sed 's/^/>/' lincRNA.genes > lincRNA.genes.modified

#Modify the TE-containing transcript list to include ">"
sed 's/^/>/' putative_intergenic.genes.not.genes.fa.blast.out > List_of_TE_containing_transcripts.txt

# Extract the sequences
python /evolinc_docker/extract_sequences-1.py lincRNA.genes.modified putative_intergenic.genes.not.genes.fa lincRNA.genes.fa

#Extract TE-containing sequences for user
python /evolinc_docker/extract_sequences-1.py List_of_TE_containing_transcripts.txt putative_intergenic.genes.not.genes.fa TE_containing_transcripts.fa
sed -i 's/_/./g' TE_containing_transcripts.fa
sed -i 's/gene=//g' TE_containing_transcripts.fa
#Create a bed file of TE-containing INTERGENIC transcripts for user
cut -f 1 -d "." putative_intergenic.genes.not.genes.fa.blast.out > TE_containing_transcript_list_transcript_ID_only.txt
grep -F -f TE_containing_transcript_list_transcript_ID_only.txt ../comparefile.gtf > TE_containing_transcripts.gtf
gff2bed < TE_containing_transcripts.gtf > TE_containing_transcripts.bed

ELAPSED_TIME_1=$(($SECONDS - $START_TIME_1))
echo "Elapsed time for step 1 is" $ELAPSED_TIME_1 "seconds" > ../$output/elapsed_time-evolinc-i.txt

# STEP 2:
START_TIME_2=$SECONDS
#Extract lincRNA candidates from original cuffmerge GTF file, using unmodified lincRNA.genes file
awk -F"." '{print $1}' lincRNA.genes > lincRNA.genes.id
sed -i 's~>~~g' lincRNA.genes.id
grep -F -f lincRNA.genes.id ../comparefile.gtf > filtered.lincRNA.gtf
gff2bed < filtered.lincRNA.gtf > lincRNA.prefilter.bed
awk 'BEGIN {OFS=FS="\t"} {gsub(/\./,"+",$6)}1' lincRNA.prefilter.bed > temp && mv temp lincRNA.prefilter.bed

ELAPSED_TIME_2=$(($SECONDS - $START_TIME_2))
echo "Elapsed time for Step 2 is" $ELAPSED_TIME_2 "seconds" >> ../$output/elapsed_time-evolinc-i.txt

# STEP 3:
START_TIME_3=$SECONDS
if [ ! -z $user_referencegff ];
then
    sed 's~^~>~g' ../$user_referencegff | sed 's~^>0*~>~g' | sed 's~^>Chr0*~>~g' | sed 's~^>Scaffold0*~>~g' | sed 's~^>~~g' > user_referencegff.gff
    intersectBed -a lincRNA.prefilter.bed -b user_referencegff.gff -u -s > SOT.genes.all.bed
else
    sed 's~^~>~g' $referencegff | sed 's~^>0*~>~g' | sed 's~^>Chr0*~>~g' | sed 's~^>Scaffold0*~>~g' | sed 's~^>~~g' > referencegff.gff
    intersectBed -a lincRNA.prefilter.bed -b referencegff.gff -u -s > SOT.genes.all.bed
fi

# Get the IDs of the overlapping exons.
cut -f 10 SOT.genes.all.bed | awk -F " " '{print $2}'| sort | uniq | sed 's~;~~g' > SOT.lincRNA.ids.txt

#create a bed file with all OTs removed
grep -vFf SOT.lincRNA.ids.txt lincRNA.prefilter.bed > lincRNA.noSOT.bed

# Create a bed file for overlapping exons from the prefilter file
grep -Ff SOT.lincRNA.ids.txt lincRNA.prefilter.bed > SOT.lincRNA.all.bed

#Create a list of all OT transcripts, including all exons related (part of the same transcript) as the exons found to be overlapping genes.
cut -f 10 SOT.lincRNA.all.bed | awk -F " " '{for(i=1;i<=NF;i++){if ($i ~/TCONS/) {print $i}}}'| sort | uniq | sed 's~;~~g' |sed 's~"~~g' | sed 's~zero_length_insertion=True~~g' |sed 's/^/>/' > SOT.lincRNA.all.txt

# Move SOT to a new file
grep -A 1 -f SOT.lincRNA.all.txt lincRNA.genes.fa > SOT.lincRNA.fa
#Clean up FASTA file
sed -i 's~--~~g' SOT.lincRNA.fa # still has an extra new line

ELAPSED_TIME_3=$(($SECONDS - $START_TIME_3))
echo "Elapsed time for Step 3 is" $ELAPSED_TIME_3 "seconds" >> ../$output/elapsed_time-evolinc-i.txt

# STEP 4:
START_TIME_4=$SECONDS
# Identify transcripts that are overlapping in the opposite direction (AOT)
if [ ! -z $user_referencegff ];
then
    sed 's~^~>~g' ../$user_referencegff | sed 's~^>0*~>~g' | sed 's~^>Chr0*~>~g' | sed 's~^>Scaffold0*~>~g' | sed 's~^>~~g' > user_referencegff.gff
    intersectBed -a lincRNA.noSOT.bed -b user_referencegff.gff -u -S > AOT.lincRNA.genes.all.bed
else
    sed 's~^~>~g' $referencegff | sed 's~^>0*~>~g' | sed 's~^>Chr0*~>~g' | sed 's~^>Scaffold0*~>~g' | sed 's~^>~~g' > referencegff.gff
    intersectBed -a lincRNA.noSOT.bed -b referencegff.gff -u -S > AOT.lincRNA.genes.all.bed
fi

# Make a list from the above file-These are the exons that overlapped
cut -f 10 AOT.lincRNA.genes.all.bed | awk -F " " '{print $2}'| sort | uniq | sed 's~;~~g' > AOT.lincRNA.ids.txt

# Generate a bed file of the entire transcript of exons that were found to be AOT
grep -Ff AOT.lincRNA.ids.txt lincRNA.noSOT.bed > AOT.lincRNA.all.bed

#Create a list of all AOT transcripts
cut -f 10 AOT.lincRNA.all.bed | awk -F " " '{for(i=1;i<=NF;i++){if ($i ~/TCONS/) {print $i}}}'| sort | uniq | sed 's~;~~g' |sed 's~"~~g' | sed 's~zero_length_insertion=True~~g' |sed 's/^/>/' > AOT.lincRNA.all.txt

#create a bed file with all AOT and SOT removed
grep -vFf AOT.lincRNA.ids.txt lincRNA.noSOT.bed > lincRNA.postfilter.bed

# Move NATs to a new file
grep -A 1 -f AOT.lincRNA.all.txt lincRNA.genes.fa > AOT.lincRNA.fa
#Clean up FASTA file
sed -i 's~--~~g' AOT.lincRNA.fa # still has an extra new line

ELAPSED_TIME_4=$(($SECONDS - $START_TIME_4))
echo "Elapsed time for Step 4 is" $ELAPSED_TIME_4 "seconds" >> ../$output/elapsed_time-evolinc-i.txt

###repeating Steps 2-4 on the overlapping transcripts to sort them into SOT or AOT

#Extract lncRNA candidates from original cuffmerge GTF file, using unmodified lncRNA.genes file
awk -F"." '{print $1}' transcripts.all.overlapping.filter.not.genes > lncRNA.genes.id
sed -i 's~>~~g' lncRNA.genes.id
grep -F -f lncRNA.genes.id ../comparefile.gtf > filtered.lncRNA.gtf
gff2bed < filtered.lncRNA.gtf > lncRNA.prefilter.bed
awk 'BEGIN {OFS=FS="\t"} {gsub(/\./,"+",$6)}1' lncRNA.prefilter.bed > temp && mv temp lncRNA.prefilter.bed

if [ ! -z $user_referencegff ];
then
    intersectBed -a lncRNA.prefilter.bed -b ../$user_referencegff -u -s > SOT.lncRNA.genes.all.bed
else   
    intersectBed -a lncRNA.prefilter.bed -b $referencegff -u -s > SOT.lncRNA.genes.all.bed
fi

# Get the IDs of the overlapping exons.
cut -f 10 SOT.lncRNA.genes.all.bed | awk -F " " '{print $2}'| sort | uniq | sed 's~;~~g' > SOT.lncRNA.ids.txt

#create a bed file with all OTs removed
grep -vFf SOT.lncRNA.ids.txt lncRNA.prefilter.bed > lncRNA.noSOT.bed

# Create a bed file for overlapping exons from the prefilter file
grep -Ff SOT.lncRNA.ids.txt lncRNA.prefilter.bed > SOT.lncRNA.all.bed

#Create a list of all OT transcripts, including all exons related (part of the same transcript) as the exons found to be overlapping genes.
cut -f 10 SOT.lncRNA.all.bed | awk -F " " '{for(i=1;i<=NF;i++){if ($i ~/TCONS/) {print $i}}}'| sort | uniq | sed 's~;~~g' |sed 's~"~~g' | sed 's~zero_length_insertion=True~~g' |sed 's/^/>/' > SOT.lncRNA.all.txt

# Move SOT to a new file
grep -A 1 -f SOT.lncRNA.all.txt transcripts.all.overlapping.filter.not.genes.fa > SOT.lncRNA.fa
#Clean up FASTA file
sed -i 's~--~~g' SOT.lncRNA.fa # still has an extra new line

# Identify transcripts that are overlapping in the opposite direction (AOT)
if [ ! -z $user_referencegff ];
then
    intersectBed -a lncRNA.noSOT.bed -b ../$user_referencegff -u -S > AOT.lncRNA.genes.all.bed
else
    intersectBed -a lncRNA.noSOT.bed -b $referencegff -u -S > AOT.lncRNA.genes.all.bed
fi

# Make a list from the above file-These are the exons that overlapped
cut -f 10 AOT.lncRNA.genes.all.bed | awk -F " " '{print $2}'| sort | uniq | sed 's~;~~g' > AOT.lncRNA.ids.txt

# Generate a bed file of the entire transcript of exons that were found to be AOT
grep -Ff AOT.lncRNA.ids.txt lncRNA.noSOT.bed > AOT.lncRNA.all.bed

#Create a list of all AOT transcripts
cut -f 10 AOT.lncRNA.all.bed | awk -F " " '{for(i=1;i<=NF;i++){if ($i ~/TCONS/) {print $i}}}'| sort | uniq | sed 's~;~~g' |sed 's~"~~g' | sed 's~zero_length_insertion=True~~g' |sed 's/^/>/' > AOT.lncRNA.all.txt

#create a bed file with all AOT and SOT removed
grep -vFf AOT.lncRNA.ids.txt lncRNA.noSOT.bed > lncRNA.postfilter.bed

# Move NATs to a new file
grep -A 1 -f AOT.lncRNA.all.txt transcripts.all.overlapping.filter.not.genes.fa > AOT.lncRNA.fa
#Clean up FASTA file
sed -i 's~--~~g' AOT.lncRNA.fa # still has an extra new line
###End of Step 2-4 repeat

#Combine SOT and AOT fasta and bed files into one set of each AOT and SOT for downstream demographics step

cat SOT.lncRNA.fa SOT.lincRNA.fa > SOT.fa
cat AOT.lncRNA.fa AOT.lincRNA.fa > AOT.fa
cat SOT.lncRNA.all.bed SOT.lincRNA.all.bed > SOT.all.bed
cat AOT.lncRNA.all.bed AOT.lincRNA.all.bed > AOT.all.bed

# STEP 5: Generating final lincRNAs.bed file and the final set of lincRNA sequences
START_TIME_5=$SECONDS
# Make a list from the lincRNA.postfilter.bed file in order to know which lincRNAs to extract from the lincRNA.genes.fa file
cut -f 10 lincRNA.postfilter.bed | awk -F " " '{for(i=1;i<=NF;i++){if ($i ~/TCONS/) {print $i ".gene=" $2}}}'| sort | uniq | sed 's~;~~g' |sed 's~"~~g' | sed 's~zero_length_insertion=True~~g' |sed 's/^/>/' > lincRNA.genes.filtered.uniq.genes
# The above line is imperfect right now, as there may be situations where the TCONS already has a gene ID that it gets called by gffread, whereas we replace with the XLOC. 
# Move lincRNA genes to a new file
python /evolinc_docker/extract_sequences-1.py lincRNA.genes.filtered.uniq.genes lincRNA.genes.fa lincRNAs.fa

# Sort the bed file for final output
sortBed -i lincRNA.postfilter.bed > lincRNAs.bed

ELAPSED_TIME_5=$(($SECONDS - $START_TIME_5))
echo "Elapsed time for Step 5 is" $ELAPSED_TIME_5 "seconds" >> ../$output/elapsed_time-evolinc-i.txt

# STEP 6: 
START_TIME_6=$SECONDS
# Update the cufflinks gtf file, only has known genes plus lincRNAs
#This first step removes all the overlapping and partial transcripts
grep -v 'class_code "x"' ../comparefile.gtf | grep -v 'class_code "s"' | grep -v 'class_code "o"' | grep -v 'class_code "e"' | grep -v 'class_code "i"' >part_modified_lincRNA.gtf
#This modifies all entries in the part_modifified.gtf file so that the second column contains lincRNA instead of Cufflinks
python /evolinc_docker/update_gtf.py lincRNAs.fa part_modified_lincRNA.gtf modified_lincRNA.updated.gtf
#Clean up quotation marks left over from update_gtf.py
sed 's~""~"~g' modified_lincRNA.updated.gtf | sed 's~"gene_id~gene_id~g' | sed 's~;"~;~g' >modified_lincRNA.updated_2.gtf
#Needed to keep "u" class codes in the annotation file until we could modify their IDs. Once we did this (above) we can now remove lines which contain BOTH cufflinks and "u"
grep -v -e 'Cufflinks.*class_code "u"' modified_lincRNA.updated_2.gtf >lincRNA.updated.gtf
rm ../comparefile.gtf
ELAPSED_TIME_6=$(($SECONDS - $START_TIME_6))
echo "Elapsed time for Step 6 is" $ELAPSED_TIME_6 "seconds" >> ../$output/elapsed_time-evolinc-i.txt

# STEP 7:
START_TIME_7=$SECONDS
# Demographics for all lincRNA
python /evolinc_docker/seq_length.py lincRNAs.fa > lincRNA.txt

cat lincRNA.txt | sort -n | cut -f 2 | awk 'NR == 1 { max=$1; min=$1; sum=0 }
   { if ($1>max) max=$1; if ($1<min) min=$1; sum+=$1;}
   END {printf "Total number of lincRNAs (including isoforms): %d\nTotal length(bp): %d\nSmallest lincRNA(bp): %d\nLargest lincRNA(bp): %d\nAverage length of lincRNA(bp): %f\n", NR, sum, min, max, sum/NR}' > lincRNA_demographics.txt

# Identify the number of unique lincRNAs in bed file, and update the demographics file
uniquelincRNAcount=$(cut -f 10 lincRNAs.bed | awk -F " " '{print $2}'| sort | uniq | grep -c "XLOC")
echo "Total number of unique lincRNAs = $uniquelincRNAcount" >> lincRNA_demographics.txt

# Demographics for SOT lncRNAs
python /evolinc_docker/seq_length.py SOT.fa > SOT.lincRNA.txt

cat SOT.lincRNA.txt | sort -n | cut -f 2 | awk 'NR == 1 { max=$1; min=$1; sum=0 }
   { if ($1>max) max=$1; if ($1<min) min=$1; sum+=$1;}
   END {printf "Total number of lincRNAs (including isoforms): %d\nTotal length(bp): %d\nSmallest lincRNA(bp): %d\nLargest lincRNA(bp): %d\nAverage length of lincRNA(bp): %f\n", NR, sum, min, max, sum/NR}' > SOT_lincRNA_demographics.txt

# Identify the number of unique lincRNAs in bed file, and update the demographics file
uniquelincRNAcount=$(cut -f 10 SOT.all.bed | awk -F " " '{print $2}'| sort | uniq | grep -c "XLOC")
echo "Total number of unique lincRNAs = $uniquelincRNAcount" >> SOT_lincRNA_demographics.txt

# # Demographics for AOT lncRNAs
python /evolinc_docker/seq_length.py AOT.fa > AOT.lincRNA.txt

cat AOT.lincRNA.txt | sort -n | cut -f 2 | awk 'NR == 1 { max=$1; min=$1; sum=0 }
   { if ($1>max) max=$1; if ($1<min) min=$1; sum+=$1;}
   END {printf "Total number of lincRNAs (including isoforms): %d\nTotal length(bp): %d\nSmallest lincRNA(bp): %d\nLargest lincRNA(bp): %d\nAverage length of lincRNA(bp): %f\n", NR, sum, min, max, sum/NR}' > AOT_lincRNA_demographics.txt

# Identify the number of unique lincRNAs in bed file, and update the demographics file
uniquelincRNAcount=$(cut -f 10 AOT.all.bed | awk -F " " '{print $2}'| sort | uniq | grep -c "XLOC")
echo "Total number of unique lincRNAs = $uniquelincRNAcount" >> AOT_lincRNA_demographics.txt

ELAPSED_TIME_7=$(($SECONDS - $START_TIME_7))
echo "Elapsed time for Step 7 is" $ELAPSED_TIME_7 "seconds" >> ../$output/elapsed_time-evolinc-i.txt

# STEP 8:
START_TIME_8=$SECONDS
# Copy the files to the outputfiles
cp lincRNAs.bed lincRNAs.fa lincRNA.updated.gtf lincRNA_demographics.txt ../$output
mkdir Other_lncRNA
cp SOT.fa SOT.all.bed AOT.fa AOT.all.bed SOT_lincRNA_demographics.txt AOT_lincRNA_demographics.txt TE_containing_transcripts.fa TE_containing_transcripts.bed Other_lncRNA
cp -r Other_lncRNA ../$output

ELAPSED_TIME_8=$(($SECONDS - $START_TIME_8))
echo "Elapsed time for Step 8 is" $ELAPSED_TIME_8 "seconds" >> ../$output/elapsed_time-evolinc-i.txt

# Optional STEP(S):
START_TIME_O1=$SECONDS

if [ ! -z $cagefile ] && [ ! -z $knownlinc ]; 
then
     sed 's~^~>~g' ../$cagefile | sed 's~^>0*~>~g' | sed 's~^>Chr0*~>~g' | sed 's~^>Scaffold0*~>~g' | sed 's~^>~~g' > cage_file &&
     gff2bed < cage_file > AnnotatedPEATPeaks.bed &&
     sortBed -i AnnotatedPEATPeaks.bed > AnnotatedPEATPeaks.sorted.bed &&
     closestBed -a lincRNAs.bed -b AnnotatedPEATPeaks.sorted.bed -s -D a > closest_output.txt && grep 'exon_number "1"' closest_output.txt > closest_output_exon_1_only.txt &&     
     python /evolinc_docker/closet_bed_compare.py closest_output_exon_1_only.txt lincRNAs.fa lincRNAs.with.CAGE.support.annotated.fa &&
     
     sed 's~^~>~g' ../$knownlinc | sed 's~^>0*~>~g' | sed 's~^>Chr0*~>~g' | sed 's~^>Scaffold0*~>~g' | sed 's~^>~~g' > knownlinc_file &&
     gff2bed < knownlinc_file > known_lncRNAs.bed &&
     sortBed -i known_lncRNAs.bed > known_lncRNAs.sorted.bed &&
     intersectBed -a lincRNAs.bed -b known_lncRNAs.sorted.bed > intersect_output.txt &&
     intersectBed -wb -a lincRNAs.bed -b known_lncRNAs.sorted.bed > intersect_output2.txt &&
     sed 's~gene_id;~gene_id ~g' intersect_output2.txt | awk -F "\t" '{print $10 ";" $20}' | sed 's~\t~~g' | awk -F ";" '{for(i=1;i<=NF;i++){if ($i ~ /gene_id/ || $i ~ /ID=/ || $i ~ /transcript_id /){print $i}}}' | sed 's~gene_id ~~g' | sed 's~"~~g' | sed 's~\n~\t~g' | sed 's~ID=~~g' | sed 's/transcript_id//' | xargs -n 3 | awk '{print $2 ".gene=" $1 "\t" $3}' | sort > temp && mv temp intersect_output2.txt
     if [ ! -s intersect_output2.txt ]; then # non-empty intersect_output file
        python /evolinc_docker/interesect_bed_compare.py intersect_output.txt lincRNAs.fa lincRNAs.overlapping.known.lincs.fa
     else # empty intersect file
        touch intersect_output.txt &&
        python /evolinc_docker/interesect_bed_compare.py intersect_output.txt lincRNAs.fa lincRNAs.overlapping.known.lincs.fa
     fi
     python /evolinc_docker/lincRNA_fig.py lincRNAs.fa lincRNAs.with.CAGE.support.annotated.fa lincRNAs.overlapping.known.lincs.fa &&
     Rscript /evolinc_docker/final_summary_table_gen_evo-I.R --lincRNA lincRNAs.fa --lincRNAbed lincRNAs.bed --overlap lincRNAs.overlapping.known.lincs.fa --tss lincRNAs.with.CAGE.support.annotated.fa &&
     sed -i 's/_/./g' lincRNAs.with.CAGE.support.annotated.fa
     sed -i 's/gene=//g' lincRNAs.with.CAGE.support.annotated.fa
     sed -i 's/_/./g' lincRNAs.overlapping.known.lincs.fa
     sed -i 's/gene=//g' lincRNAs.overlapping.known.lincs.fa
     cp lincRNAs.with.CAGE.support.annotated.fa lincRNAs.overlapping.known.lincs.fa lincRNA_piechart.png final_Summary_table.tsv ../$output

elif [ ! -z $cagefile ]; 
then
     sed 's~^~>~g' ../$cagefile | sed 's~^>0*~>~g' | sed 's~^>Chr0*~>~g' | sed 's~^>Scaffold0*~>~g' | sed 's~^>~~g' > cage_file &&
     gff2bed < cage_file > AnnotatedPEATPeaks.bed &&
     sortBed -i AnnotatedPEATPeaks.bed > AnnotatedPEATPeaks.sorted.bed &&
     closestBed -a lincRNAs.bed -b AnnotatedPEATPeaks.sorted.bed -s -D a > closest_output.txt && grep 'exon_number "1"' closest_output.txt > closest_output_exon_1_only.txt &&     
     python /evolinc_docker/closet_bed_compare.py closest_output_exon_1_only.txt lincRNAs.fa lincRNAs.with.CAGE.support.annotated.fa &&
     Rscript /evolinc_docker/final_summary_table_gen_evo-I.R --lincRNA lincRNAs.fa --lincRNAbed lincRNAs.bed --tss lincRNAs.with.CAGE.support.annotated.fa &&
     sed -i 's/_/./g' lincRNAs.with.CAGE.support.annotated.fa
     sed -i 's/gene=//g' lincRNAs.with.CAGE.support.annotated.fa
     cp lincRNAs.with.CAGE.support.annotated.fa final_Summary_table.tsv ../$output

elif [ ! -z $knownlinc ]; 
then
     sed 's~^~>~g' ../$knownlinc | sed 's~^>0*~>~g' | sed 's~^>Chr0*~>~g' | sed 's~^>Scaffold0*~>~g' | sed 's~^>~~g' > knownlinc_file &&
     gff2bed < knownlinc_file > known_lncRNAs.bed &&
     sortBed -i known_lncRNAs.bed > known_lncRNAs.sorted.bed &&
     intersectBed -a lincRNAs.bed -b known_lncRNAs.sorted.bed > intersect_output.txt &&
     intersectBed -wb -a lincRNAs.bed -b known_lncRNAs.sorted.bed > intersect_output2.txt &&
     sed 's~gene_id;~gene_id ~g' intersect_output2.txt | awk -F "\t" '{print $10 ";" $20}' | sed 's~\t~~g' | awk -F ";" '{for(i=1;i<=NF;i++){if ($i ~ /gene_id/ || $i ~ /ID=/ || $i ~ /transcript_id /){print $i}}}' | sed 's~gene_id ~~g' | sed 's~"~~g' | sed 's~\n~\t~g' | sed 's~ID=~~g' | sed 's/transcript_id//' | xargs -n 3 | awk '{print $2 ".gene=" $1 "\t" $3}' | sort > temp && mv temp intersect_output2.txt
     if [ ! -s intersect_output2.txt ]; then # non-empty intersect_output file
        python /evolinc_docker/interesect_bed_compare.py intersect_output.txt lincRNAs.fa lincRNAs.overlapping.known.lincs.fa
     else # empty intersect file
        touch intersect_output.txt &&
        python /evolinc_docker/interesect_bed_compare.py intersect_output.txt lincRNAs.fa lincRNAs.overlapping.known.lincs.fa
     fi
     Rscript /evolinc_docker/final_summary_table_gen_evo-I.R --lincRNA lincRNAs.fa --lincRNAbed lincRNAs.bed --overlap lincRNAs.overlapping.known.lincs.fa &&
     sed -i 's/_/./g' lincRNAs.overlapping.known.lincs.fa
     sed -i 's/gene=//g' lincRNAs.overlapping.known.lincs.fa
     cp lincRNAs.overlapping.known.lincs.fa final_Summary_table.tsv ../$output

else
    Rscript /evolinc_docker/final_summary_table_gen_evo-I.R --lincRNA lincRNAs.fa --lincRNAbed lincRNAs.bed
    cp final_Summary_table.tsv ../$output
fi

ELAPSED_TIME_O1=$(($SECONDS - $START_TIME_O1))
echo "Elapsed time for Optional Step(s) is" $ELAPSED_TIME_O1 "seconds" >> ../$output/elapsed_time-evolinc-i.txt

# # remove all the other files
rm -r ../transcripts_u_filter.fa.transdecoder_dir
rm ../*.fa*.*

### Clean up fasta headers
cd ../$output
sed -i 's/_/./g' lincRNAs.fa
sed -i 's/gene=//g' lincRNAs.fa
sed -i 's/_/./g' Other_lncRNA/AOT.fa
sed -i 's/gene=//g' Other_lncRNA/AOT.fa
sed -i 's/_/./g' Other_lncRNA/SOT.fa
sed -i 's/gene=//g' Other_lncRNA/SOT.fa
echo "All necessary files written to" $output
echo "Finished Evolinc-part-I!"
echo `date`

ELAPSED_TIME=$(($SECONDS - $START_TIME))
echo "Total elapsed time is" $ELAPSED_TIME "seconds" >> ../$output/elapsed_time-evolinc-i.txt
