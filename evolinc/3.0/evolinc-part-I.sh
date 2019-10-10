#!/bin/bash
# Upendra Kumar Devisetty
# Script to process cuffcompare output file to generate long non-coding RNA
# Usage: 
# sh evolinc-part-I.sh -c sample.data/cuffcompare_out_annot_no_annot.combined.gtf -g sample.data/Brassica_rapa_v1.2_genome.fa -r sample.data/Brassica_rapa_v1.2.gff -b sample.data/TE_RNA_transcripts.fa -t CAGE_file.txt -x 2> errors.txt > output.txt

usage() {
      echo ""
      echo "Usage : sh $0 -c cuffcompare -g genome -r gff -o output -threads [-b TE_RNA] [-t CAGE_RNA] [-x Known_lincRNA]"
      echo ""

cat <<'EOF'

  -c </path/to/cuffcompare output file>

  -g </path/to/reference genome file>

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

while getopts ":b:c:g:hr:t:x:o:n:" opt; do
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

# Creat a directory to move all the output files
mkdir $output


# STEP 1:
START_TIME_1=$SECONDS
# Extracting classcode u transcripts, making fasta file, removing transcripts > 200 and selecting protein coding transcripts
grep -e '"u"' -e '"x"' -e '"s"' -e '"o"' -e '"e"' -e '"i"' $comparefile | gffread -w \
    transcripts_u.fa -g $referencegenome - && python /evolinc_docker/get_gene_length_filter.py transcripts_u.fa \
    transcripts_u_filter.fa

# Make directory
mkdir transcripts_u_filter.fa.transdecoder_dir

# For files that are too large, this script will split them up into 1000 sequences each  
perl /evolinc_docker/split_multifasta.pl --input_file transcripts_u_filter.fa --output_dir transcripts_u_filter.fa.transdecoder_dir --seqs_per_file=1000

# Modifying the header
sed 's/ .*//' transcripts_u_filter.fa | sed -ne 's/>//p' > transcripts_u_filter.fa.genes

# Move the transcript files to this directory transcripts_u_filter.fa.transdecoder_dir
mv transcripts_u.fa transcripts_u_filter.fa transcripts_u_filter.fa.genes transcripts_u_filter.fa.transdecoder_dir/

# Change the directory
cd transcripts_u_filter.fa.transdecoder_dir

# Remove empty lines and change the file suffix
for i in *.fsa; do sed '/^\s*$/d' $i > $i.fasta; done
rm *.fsa
for i in *fasta; do mv $i "`basename $i .fsa.fasta`.fasta"; done

# Run transdecoder now
for files in *fasta; do TransDecoder.LongOrfs -t $files; done

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
grep -v -F -f longest_orfs_cat.cds.pep.blastp.genes transcripts_u_filter.fa.genes > transcripts_u_filter.not.genes #I added the -F here, it speeds things up quite a bit as it is searching for exact strings.

sed 's/^/>/' transcripts_u_filter.not.genes > temp && mv temp transcripts_u_filter.not.genes # changed name here 

# Extract fasta file
python /evolinc_docker/extract_sequences.py transcripts_u_filter.not.genes transcripts_u_filter.fa transcripts_u_filter.not.genes.fa 

sed 's/ /./' transcripts_u_filter.not.genes.fa > temp && mv temp transcripts_u_filter.not.genes.fa

# Blast the fasta file to TE RNA db
if [ ! -z $blastfile ]; then
     makeblastdb -in ../$blastfile -dbtype nucl -out ../$blastfile.blast.out &&
     blastn -query transcripts_u_filter.not.genes.fa -db ../$blastfile.blast.out -out transcripts_u_filter.not.genes.fa.blast.out -outfmt 6 -num_threads $threads # no blast hits here
else
    touch transcripts_u_filter.not.genes.fa.blast.out
fi

# Filter the output to select the best transcript based on e-value and bit-score
python /evolinc_docker/filter_sequences.py transcripts_u_filter.not.genes.fa.blast.out transcripts_u_filter.not.genes.fa.blast.out.filtered

# Modify the header in the fasta file to extract header only
grep ">" transcripts_u_filter.not.genes.fa | sed 's/>//' > transcripts_u_filter.not.genes_only

# Now remove the blast hits from the fasta file
python /evolinc_docker/fasta_remove.py transcripts_u_filter.not.genes.fa.blast.out.filtered transcripts_u_filter.not.genes_only lincRNA.genes

# Modify the fasta header to include ">", generating a new file so as to keep the lincRNA.genes file intact for later use.
sed 's/^/>/' lincRNA.genes > lincRNA.genes.modified

#Modify the TE-containing transcript list to include ">"
sed 's/^/>/' transcripts_u_filter.not.genes.fa.blast.out.filtered > List_of_TE_containing_transcripts.txt

# Extract the sequences
python /evolinc_docker/extract_sequences-1.py lincRNA.genes.modified transcripts_u_filter.not.genes.fa lincRNA.genes.fa

#Extract TE-containing sequences for user
python /evolinc_docker/extract_sequences-1.py List_of_TE_containing_transcripts.txt transcripts_u_filter.not.genes.fa TE_containing_transcripts.fa

#Create a bed file of TE-containing transcripts for user
cut -f 1 -d "." transcripts_u_filter.not.genes.fa.blast.out.filtered > TE_containing_transcript_list_transcript_ID_only.txt
grep -F -f TE_containing_transcript_list_transcript_ID_only.txt ../$comparefile > TE_containing_transcripts.gtf
gff2bed < TE_containing_transcripts.gtf > TE_containing_transcripts.bed

ELAPSED_TIME_1=$(($SECONDS - $START_TIME_1))
echo "Elapsed time for step 1 is" $ELAPSED_TIME_1 "seconds" > ../$output/elapsed_time-evolinc-i.txt

# STEP 2:
START_TIME_2=$SECONDS
#Extract lincRNA candidates from original cuffmerge GTF file, using unmodified lincRNA.genes file
awk -F"." '{print $1}' lincRNA.genes > lincRNA.genes.id
grep -F -f lincRNA.genes.id ../$comparefile > filtered.lincRNA.gtf
gff2bed < filtered.lincRNA.gtf > lincRNA.prefilter.bed
awk 'BEGIN {OFS=FS="\t"} {gsub(/\./,"+",$6)}1' lincRNA.prefilter.bed > temp && mv temp lincRNA.prefilter.bed

ELAPSED_TIME_2=$(($SECONDS - $START_TIME_2))
echo "Elapsed time for Step 2 is" $ELAPSED_TIME_2 "seconds" >> ../$output/elapsed_time-evolinc-i.txt

# STEP 3:
START_TIME_3=$SECONDS
intersectBed -a lincRNA.prefilter.bed -b ../$referencegff -u -s > SOT.genes.all.bed

# Get the IDs of the overlapping exons.
cut -f 10 SOT.genes.all.bed | awk -F " " '{print $2}'| sort | uniq | sed 's~;~~g' > SOT.ids.txt

#create a bed file with all OTs removed
grep -vFf SOT.ids.txt lincRNA.prefilter.bed > lincRNA.noSOT.bed

# Create a bed file for overlapping exons from the prefilter file
grep -Ff SOT.ids.txt lincRNA.prefilter.bed > SOT.all.bed

#Create a list of all OT transcripts, including all exons related (part of the same transcript) as the exons found to be overlapping genes.
cut -f 10 SOT.all.bed | awk -F " " '{for(i=1;i<=NF;i++){if ($i ~/TCONS/) {print $i}}}'| sort | uniq | sed 's~;~~g' |sed 's~"~~g' | sed 's~zero_length_insertion=True~~g' |sed 's/^/>/' > SOT.all.txt

# Move SOT to a new file
#python /evolinc_docker/extract_sequences-1.py SOT.all.txt lincRNA.genes.fa SOT.fa
grep -A 1 -f SOT.all.txt lincRNA.genes.fa > SOT.fa
#Clean up FASTA file
sed -i 's~--~~g' SOT.fa # still has an extra new line
ELAPSED_TIME_3=$(($SECONDS - $START_TIME_3))
echo "Elapsed time for Step 3 is" $ELAPSED_TIME_3 "seconds" >> ../$output/elapsed_time-evolinc-i.txt

# STEP 4:
START_TIME_4=$SECONDS
# Identify transcripts that are overlapping in the opposite direction (AOT)
intersectBed -a lincRNA.noSOT.bed -b ../$referencegff -u -S > AOT.genes.all.bed

# Make a list from the above file-These are the exons that overlapped
cut -f 10 AOT.genes.all.bed | awk -F " " '{print $2}'| sort | uniq | sed 's~;~~g' > AOT.ids.txt

# Generate a bed file of the entire transcript of exons that were found to be AOT
grep -Ff AOT.ids.txt lincRNA.noSOT.bed > AOT.all.bed

#Create a list of all AOT transcripts
cut -f 10 AOT.all.bed | awk -F " " '{for(i=1;i<=NF;i++){if ($i ~/TCONS/) {print $i}}}'| sort | uniq | sed 's~;~~g' |sed 's~"~~g' | sed 's~zero_length_insertion=True~~g' |sed 's/^/>/' > AOT.all.txt

#create a bed file with all AOT and SOT removed
grep -vFf AOT.ids.txt lincRNA.noSOT.bed > lincRNA.postfilter.bed

# Move NATs to a new file
#python /evolinc_docker/extract_sequences-1.py AOT.all.txt lincRNA.genes.fa AOT.fa
grep -A 1 -f AOT.all.txt lincRNA.genes.fa > AOT.fa
#Clean up FASTA file
sed -i 's~--~~g' AOT.fa # still has an extra new line

ELAPSED_TIME_4=$(($SECONDS - $START_TIME_4))
echo "Elapsed time for Step 4 is" $ELAPSED_TIME_4 "seconds" >> ../$output/elapsed_time-evolinc-i.txt

# STEP 5: Generating final lincRNA.bed file and the final set of lincRNA sequences
START_TIME_5=$SECONDS
# Make a list from the lincRNA.postfilter.bed file in order to know which lincRNAs to extract from the lincRNA.genes.fa file
cut -f 10 lincRNA.postfilter.bed | awk -F " " '{for(i=1;i<=NF;i++){if ($i ~/TCONS/) {print $i ".gene=" $2}}}'| sort | uniq | sed 's~;~~g' |sed 's~"~~g' | sed 's~zero_length_insertion=True~~g' |sed 's/^/>/' > lincRNA.genes.filtered.uniq.genes
# The above line is imperfect right now, as there may be situations where the TCONS already has a gene ID that it gets called by gffread, whereas we replace with the XLOC. 
# Move lincRNA genes to a new file
python /evolinc_docker/extract_sequences-1.py lincRNA.genes.filtered.uniq.genes lincRNA.genes.fa All.lincRNAs.fa

# Sort the bed file for final output
sortBed -i lincRNA.postfilter.bed > lincRNA.bed

ELAPSED_TIME_5=$(($SECONDS - $START_TIME_5))
echo "Elapsed time for Step 5 is" $ELAPSED_TIME_5 "seconds" >> ../$output/elapsed_time-evolinc-i.txt

# STEP 6: 
START_TIME_6=$SECONDS
# Update the cufflinks gtf file
python /evolinc_docker/update_gtf.py All.lincRNAs.fa ../$comparefile lincRNA.updated.gtf

ELAPSED_TIME_6=$(($SECONDS - $START_TIME_6))
echo "Elapsed time for Step 6 is" $ELAPSED_TIME_6 "seconds" >> ../$output/elapsed_time-evolinc-i.txt

# STEP 7:
START_TIME_7=$SECONDS
# Demographics for all lincRNA
quast.py All.lincRNAs.fa -R ../$referencegenome -G ../$referencegff --threads $threads -o lincRNA_demographics
sed 's/contig/lincRNA/g' lincRNA_demographics/report.txt > temp && mv temp lincRNA_demographics/report.txt
sed 's/contig/lincRNA/g' lincRNA_demographics/icarus.html > temp && mv temp lincRNA_demographics/icarus.html

# Identify the number of unique lincRNAs in bed file, and update the demographics file
uniquelincRNAcount=$(cut -f 10 lincRNA.bed | awk -F " " '{print $2}'| sort | uniq | grep -c "XLOC")
sed -i "4i Unique lincRNAs          $uniquelincRNAcount" lincRNA_demographics/report.txt
sed -i "s~# lincRNAs (>~# of total lincRNAs (including isoforms) (>~g" lincRNA_demographics/report.txt

#Clip some of the lincRNA_demographics report info and move to output folder
head -n 22 lincRNA_demographics/report.txt > lincRNA_demographics.txt
grep -v "Assembly" lincRNA_demographics.txt > temp && mv temp lincRNA_demographics.txt
cp lincRNA_demographics.txt ../$output

# Demographics for SOT lncRNAs
quast.py SOT.fa -R ../$referencegenome -G ../$referencegff --threads $threads -o SOT_demographics
sed 's/contig/lincRNA/g' SOT_demographics/report.txt > temp && mv temp SOT_demographics/report.txt
sed 's/contig/SOT/g' SOT_demographics/icarus.html > temp && mv temp SOT_demographics/icarus.html

# Modify demographics file to include number of unique transcripts (unique XLOCs)
uniqueOTcount=$(cut -f 10 SOT.all.bed | awk -F " " '{print $2}'| sort | uniq | grep -c "XLOC")
sed -i "4i Unique SOT lncRNAs           $uniqueOTcount" SOT_demographics/report.txt

# Change the terminology to reflect SOT
sed -i "s~# lincRNAs (>~# of total SOT lncRNAs (including isoforms) (>~g" SOT_demographics/report.txt
sed -i "s~lincRNA~SOT lncRNA~g" SOT_demographics/report.txt

#Clip some of the SOT_demographics report info and move to output folder
head -n 22 SOT_demographics/report.txt > SOT_demographics.txt
grep -v "Assembly" SOT_demographics.txt > temp && mv temp SOT_demographics.txt
cp SOT_demographics.txt ../$output

# Demographics for AOT lncRNAs
quast.py AOT.fa -R ../$referencegenome -G ../$referencegff --threads $threads -o AOT_demographics
sed -i 's/contig/lincRNA/g' AOT_demographics/report.txt > AOT_demographics/report.txt
sed -i 's/contig/AOT/g' AOT_demographics/icarus.html > AOT_demographics/icarus.html

# Modify demographics file to include number of unique transcripts (unique XLOCs)
uniqueNATcount=$(cut -f 10 AOT.all.bed | awk -F " " '{print $2}'| sort | uniq | grep -c "XLOC") #What if the gtf file they input doesn't use the XLOC nomenclature?
sed -i "4i Unique AOT lncRNAs        $uniqueNATcount" AOT_demographics/report.txt

# Change the terminology to reflect AOT
sed -i "s~# lincRNAs (>~# of total AOT lncRNAs (including isoforms) (>~g" AOT_demographics/report.txt
sed -i "s~lincRNA~AOT lncRNA~g" AOT_demographics/report.txt

#Clip some of the AOT_demographics report info and move to output folder
head -n 22 AOT_demographics/report.txt > AOT_demographics.txt
grep -v "Assembly" AOT_demographics.txt > temp && mv temp AOT_demographics.txt
cp AOT_demographics.txt ../$output

ELAPSED_TIME_7=$(($SECONDS - $START_TIME_7))
echo "Elapsed time for Step 7 is" $ELAPSED_TIME_7 "seconds" >> ../$output/elapsed_time-evolinc-i.txt

# STEP 8:
START_TIME_8=$SECONDS
# Copy the files to the outputfiles
cp -r lincRNA.bed All.lincRNAs.fa lincRNA_demographics lincRNA.updated.gtf SOT.fa SOT.all.bed SOT_demographics AOT.fa AOT.all.bed AOT_demographics TE_containing_transcripts.fa TE_containing_transcripts.bed ../$output

ELAPSED_TIME_8=$(($SECONDS - $START_TIME_8))
echo "Elapsed time for Step 8 is" $ELAPSED_TIME_8 "seconds" >> ../$output/elapsed_time-evolinc-i.txt

# Optional STEP(S):
START_TIME_O1=$SECONDS

if [ ! -z $cagefile ] && [ ! -z $knownlinc ]; 
then
      
     gff2bed < ../$cagefile > AnnotatedPEATPeaks.bed &&
     sortBed -i AnnotatedPEATPeaks.bed > AnnotatedPEATPeaks.sorted.bed &&
     closestBed -a lincRNA.bed -b AnnotatedPEATPeaks.sorted.bed -s -D a > closest_output.txt && grep 'exon_number "1"' closest_output.txt > closest_output_exon_1_only.txt &&     
     python /evolinc_docker/closet_bed_compare.py closest_output_exon_1_only.txt All.lincRNAs.fa lincRNAs.with.CAGE.support.annotated.fa &&
     
     gff2bed < ../$knownlinc > known_lncRNAs.bed &&
     sortBed -i known_lncRNAs.bed > known_lncRNAs.sorted.bed &&
     intersectBed -a lincRNA.bed -b known_lncRNAs.sorted.bed > intersect_output.txt &&
     intersectBed -wb -a lincRNA.bed -b known_lncRNAs.sorted.bed > intersect_output2.txt &&
     sed 's~gene_id;~gene_id ~g' intersect_output2.txt | awk -F "\t" '{print $10 ";" $20}' | sed 's~\t~~g' | awk -F ";" '{for(i=1;i<=NF;i++){if ($i ~ /gene_id/ || $i ~ /ID=/ || $i ~ /transcript_id /){print $i}}}' | sed 's~gene_id ~~g' | sed 's~"~~g' | sed 's~\n~\t~g' | sed 's~ID=~~g' | sed 's/transcript_id//' | xargs -n 3 | awk '{print $2 ".gene=" $1 "\t" $3}' | sort > temp && mv temp intersect_output2.txt
     if [ ! -s intersect_output2.txt ]; then # non-empty intersect_output file
        python /evolinc_docker/interesect_bed_compare.py intersect_output.txt All.lincRNAs.fa lincRNAs.overlapping.known.lincs.fa
     else # empty intersect file
        touch intersect_output.txt &&
        python /evolinc_docker/interesect_bed_compare.py intersect_output.txt All.lincRNAs.fa lincRNAs.overlapping.known.lincs.fa
     fi
     python /evolinc_docker/lincRNA_fig.py All.lincRNAs.fa lincRNAs.with.CAGE.support.annotated.fa lincRNAs.overlapping.known.lincs.fa &&
     Rscript /evolinc_docker/final_summary_table_gen_evo-I.R --lincRNA All.lincRNAs.fa --lincRNAbed lincRNA.bed --overlap lincRNAs.overlapping.known.lincs.fa --tss lincRNAs.with.CAGE.support.annotated.fa &&
     
     cp lincRNAs.with.CAGE.support.annotated.fa lincRNAs.overlapping.known.lincs.fa lincRNA_piechart.png final_Summary_table.tsv ../$output

elif [ ! -z $cagefile ]; 
then
     gff2bed < ../$cagefile > AnnotatedPEATPeaks.bed &&
     sortBed -i AnnotatedPEATPeaks.bed > AnnotatedPEATPeaks.sorted.bed &&
     closestBed -a lincRNA.bed -b AnnotatedPEATPeaks.sorted.bed -s -D a > closest_output.txt && grep 'exon_number "1"' closest_output.txt > closest_output_exon_1_only.txt &&     
     python /evolinc_docker/closet_bed_compare.py closest_output_exon_1_only.txt All.lincRNAs.fa lincRNAs.with.CAGE.support.annotated.fa &&
     Rscript /evolinc_docker/final_summary_table_gen_evo-I.R --lincRNA All.lincRNAs.fa --lincRNAbed lincRNA.bed --tss lincRNAs.with.CAGE.support.annotated.fa &&
     cp lincRNAs.with.CAGE.support.annotated.fa final_Summary_table.tsv ../$output

elif [ ! -z $knownlinc ]; 
then
    gff2bed < ../$knownlinc > known_lncRNAs.bed &&
     sortBed -i known_lncRNAs.bed > known_lncRNAs.sorted.bed &&
     intersectBed -a lincRNA.bed -b known_lncRNAs.sorted.bed > intersect_output.txt &&
     intersectBed -wb -a lincRNA.bed -b known_lncRNAs.sorted.bed > intersect_output2.txt &&
     sed 's~gene_id;~gene_id ~g' intersect_output2.txt | awk -F "\t" '{print $10 ";" $20}' | sed 's~\t~~g' | awk -F ";" '{for(i=1;i<=NF;i++){if ($i ~ /gene_id/ || $i ~ /ID=/ || $i ~ /transcript_id /){print $i}}}' | sed 's~gene_id ~~g' | sed 's~"~~g' | sed 's~\n~\t~g' | sed 's~ID=~~g' | sed 's/transcript_id//' | xargs -n 3 | awk '{print $2 ".gene=" $1 "\t" $3}' | sort > temp && mv temp intersect_output2.txt
     if [ ! -s intersect_output2.txt ]; then # non-empty intersect_output file
        python /evolinc_docker/interesect_bed_compare.py intersect_output.txt All.lincRNAs.fa lincRNAs.overlapping.known.lincs.fa
     else # empty intersect file
        touch intersect_output.txt &&
        python /evolinc_docker/interesect_bed_compare.py intersect_output.txt All.lincRNAs.fa lincRNAs.overlapping.known.lincs.fa
     fi
     Rscript /evolinc_docker/final_summary_table_gen_evo-I.R --lincRNA All.lincRNAs.fa --lincRNAbed lincRNA.bed --overlap lincRNAs.overlapping.known.lincs.fa &&
     cp lincRNAs.overlapping.known.lincs.fa final_Summary_table.tsv ../$output

fi

ELAPSED_TIME_O1=$(($SECONDS - $START_TIME_O1))
echo "Elapsed time for Optional Step(s) is" $ELAPSED_TIME_O1 "seconds" >> ../$output/elapsed_time-evolinc-i.txt

# # Optional STEP - 1:
# START_TIME_O1=$SECONDS
# # CAGE data
# if [ ! -z $cagefile ]; then
#      gff2bed < ../$cagefile > AnnotatedPEATPeaks.bed &&
#      sortBed -i AnnotatedPEATPeaks.bed > AnnotatedPEATPeaks.sorted.bed &&
#      closestBed -a lincRNA.bed -b AnnotatedPEATPeaks.sorted.bed -s -D a > closest_output.txt && grep 'exon_number "1"' closest_output.txt > closest_output_exon_1_only.txt &&     
#      python /evolinc_docker/closet_bed_compare.py closest_output_exon_1_only.txt All.lincRNAs.fa lincRNAs.with.CAGE.support.annotated.fa &&
#      Rscript /evolinc_docker/final_summary_table_gen_evo-I.R --lincRNA All.lincRNAs.fa --lincRNAbed lincRNA.bed --tss lincRNAs.with.CAGE.support.annotated.fa &&
#      cp lincRNAs.with.CAGE.support.annotated.fa final_Summary_table.tsv ../$output
# fi

# ELAPSED_TIME_O1=$(($SECONDS - $START_TIME_O1))
# echo "Elapsed time for Optional Step 1 is" $ELAPSED_TIME_O1 "seconds" >> ../$output/elapsed_time-evolinc-i.txt

# # Optional STEP - 2:
# START_TIME_O2=$SECONDS
# # Known lincRNA
# if [ ! -z $knownlinc ]; then
#      gff2bed < ../$knownlinc > known_lncRNAs.bed &&
#      sortBed -i known_lncRNAs.bed > known_lncRNAs.sorted.bed &&
#      intersectBed -wb -a lincRNA.bed -b known_lncRNAs.sorted.bed > intersect_output.txt &&
#      sed 's~gene_id;~gene_id ~g' intersect_output.txt | awk -F "\t" '{print $10 ";" $20}' | sed 's~\t~~g' | awk -F ";" '{for(i=1;i<=NF;i++){if ($i ~ /gene_id/ || $i ~ /ID=/ || $i ~ /transcript_id /){print $i}}}' | sed 's~gene_id ~~g' | sed 's~"~~g' | sed 's~\n~\t~g' | sed 's~ID=~~g' | sed 's/transcript_id//' | xargs -n 3 | awk '{print $2 ".gene=" $1 "\t" $3}' | sort > temp && mv temp intersect_output.txt
#      if [ ! -s intersect_output.txt ]; then # non-empty intersect_output file
#         python /evolinc_docker/interesect_bed_compare.py intersect_output.txt All.lincRNAs.fa lincRNAs.overlapping.known.lincs.fa &&
#         Rscript /evolinc_docker/final_summary_table_gen_evo-I.R --lincRNA All.lincRNAs.fa --lincRNAbed lincRNA.bed --overlap lincRNAs.overlapping.known.lincs.fa &&
#         cp lincRNAs.overlapping.known.lincs.fa final_Summary_table.tsv ../$output
#      else # empty intersect file
#         touch intersect_output.txt
#         python /evolinc_docker/interesect_bed_compare.py intersect_output.txt All.lincRNAs.fa lincRNAs.overlapping.known.lincs.fa &&
#         cp lincRNAs.overlapping.known.lincs.fa ../$output
#      fi
          
# fi

# ELAPSED_TIME_O2=$(($SECONDS - $START_TIME_O2))
# echo "Elapsed time for Optional Step 2 is" $ELAPSED_TIME_O2 "seconds" >> ../$output/elapsed_time-evolinc-i.txt

# # Pie chart if both CAGE and Knownlinc are given
# if [ ! -z $cagefile ] && [ ! -z $knownlinc ] ; then
#    python /evolinc_docker/lincRNA_fig.py All.lincRNAs.fa lincRNAs.with.CAGE.support.annotated.fa lincRNAs.overlapping.known.lincs.fa &&
#    Rscript /evolinc_docker/final_summary_table_gen_evo-I.R --lincRNA All.lincRNAs.fa --lincRNAbed lincRNA.bed --overlap lincRNAs.overlapping.known.lincs.fa --tss lincRNAs.with.CAGE.support.annotated.fa
#    cp lincRNA_piechart.png final_Summary_table.tsv ../$output
# fi

# remove all the other files
#rm -r ../transcripts_u_filter.fa.transdecoder_dir
#rm ../*.fa*.*

echo "All necessary files written to" $output
echo "Finished Evolinc-part-I!"
echo `date`

ELAPSED_TIME=$(($SECONDS - $START_TIME))
echo "Total elapsed time is" $ELAPSED_TIME "seconds" >> ../$output/elapsed_time-evolinc-i.txt