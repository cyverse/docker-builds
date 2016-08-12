#!/bin/bash
# Upendra Kumar Devisetty
# Script to process cuffcompare output file to generate long non-coding RNA
# Usage: 
# sh evolinc-part-I.sh -c sample.data/cuffcompare_out_annot_no_annot.combined.gtf -g sample.data/Brassica_rapa_v1.2_genome.fa -r sample.data/Brassica_rapa_v1.2.gff -b sample.data/TE_RNA_transcripts.fa -t CAGE_file.txt -x 2> errors.txt > output.txt

usage() {
      echo ""
      echo "Usage : sh $0 -c cuffcompare -g genome -r gff -o output [-b TE_RNA] [-t CAGE_RNA] [-x Known_lincRNA]"
      echo ""

cat <<'EOF'
  -c </path/to/cuffcompare output file>

  -g </path/to/reference genome file>

  -r </path/to/reference annotation file>

  -b </path/to/Transposable Elements file>

  -o </path/to/output file>

  -t </path/to/CAGE RNA file>
  
  -x </path/to/Known lincRNA file>

  -h Show this usage information

EOF
    exit 0
}

while getopts ":b:c:g:hr:t:x:o:" opt; do
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

# For files that are too large, this script will split them up into 100 sequences each  
perl /evolinc_docker/split_multifasta.pl --input_file transcripts_u_filter.fa --output_dir transcripts_u_filter.fa.transdecoder_dir --seqs_per_file=100

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

#This groups all the longest_orfs.cds files into one, in the transdecoder file.
find . -type f -name longest_orfs.cds -exec cat '{}' \; | cat > longest_orfs_cat.cds 

# Genes in the protein coding genes
sed 's/|.*//' longest_orfs_cat.cds | sed -ne 's/>//p' | uniq > longest_orfs.cds.genes

# Remove these protein coding genes from the filter file
grep -v -F -f longest_orfs.cds.genes transcripts_u_filter.fa.genes > transcripts_u_filter.not.genes #I added the -F here, it speeds things up quite a bit as it is searching for exact strings.

sed 's/^/>/' transcripts_u_filter.not.genes > temp && mv temp transcripts_u_filter.not.genes # changed name here 

# Extract fasta file
python /evolinc_docker/extract_sequences.py transcripts_u_filter.not.genes transcripts_u_filter.fa transcripts_u_filter.not.genes.fa 
sed 's/ /./' transcripts_u_filter.not.genes.fa > temp && mv temp transcripts_u_filter.not.genes.fa

# Blast the fasta file to TE RNA db
if [ ! -z $blastfile ]; then
     makeblastdb -in ../$blastfile -dbtype nucl -out ../$blastfile.blast.out &&
     blastn -query transcripts_u_filter.not.genes.fa -db ../$blastfile.blast.out -out transcripts_u_filter.not.genes.fa.blast.out -outfmt 6 # no blast hits her
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

# Extract the sequences
python /evolinc_docker/extract_sequences-1.py lincRNA.genes.modified transcripts_u_filter.not.genes.fa lincRNA.genes.fa

ELAPSED_TIME_1=$(($SECONDS - $START_TIME_1))
echo "Elapsed time for step 1 is" $ELAPSED_TIME_1 "seconds" > ../$output/elapsed_time-evolinc-i.txt

# STEP 2:
START_TIME_2=$SECONDS
#Extract lincRNA candidates from original cuffmerge GTF file, using unmodified lincRNA.genes file
awk -F"." '{print $1}' lincRNA.genes >lincRNA.genes.id
grep -F -f lincRNA.genes.id ../$comparefile > filtered.lincRNA.gtf
gff2bed < filtered.lincRNA.gtf > lincRNA.prefilter.bed
awk 'BEGIN {OFS=FS="\t"} {gsub(/\./,"+",$6)}1' lincRNA.prefilter.bed > temp && mv temp lincRNA.prefilter.bed

ELAPSED_TIME_2=$(($SECONDS - $START_TIME_2))
echo "Elapsed time for Step 2 is" $ELAPSED_TIME_2 "seconds" >> ../$output/elapsed_time-evolinc-i.txt

# STEP 3:
START_TIME_3=$SECONDS
intersectBed -a lincRNA.prefilter.bed -b ../$referencegff -u -s > OT.genes.all.bed

# Get the IDs of the overlapping exons.
cut -f 10 OT.genes.all.bed | awk -F " " '{print $2}'| sort | uniq | sed 's~;~~g' > OT.ids.txt

#create a bed file with all OTs removed
grep -vFf OT.ids.txt lincRNA.prefilter.bed > lincRNA.noOTs.bed

# Create a bed file for overlapping exons from the prefilter file
grep -Ff OT.ids.txt lincRNA.prefilter.bed > OT.all.bed

#Create a list of all OT transcripts, including all exons related (part of the same transcript) as the exons found to be overlapping genes.
cut -f 10 OT.all.bed | awk -F " " '{for(i=1;i<=NF;i++){if ($i ~/TCONS/) {print $i ".gene=" $2}}}'| sort | uniq | sed 's~;~~g' |sed 's~"~~g' | sed 's~zero_length_insertion=True~~g' |sed 's/^/>/' > OT.all.txt

# Move OTs to a new file
python /evolinc_docker/extract_sequences-1.py OT.all.txt lincRNA.genes.fa Overlapping.transcripts.fa

ELAPSED_TIME_3=$(($SECONDS - $START_TIME_3))
echo "Elapsed time for Step 3 is" $ELAPSED_TIME_3 "seconds" >> ../$output/elapsed_time-evolinc-i.txt

# STEP 4:
START_TIME_4=$SECONDS
# Identify transcripts that are overlapping in the opposite direction (NAT)
intersectBed -a lincRNA.noOTs.bed -b ../$referencegff -u -S > NAT.genes.all.bed

# Make a list from the above file-These are the exons that overlapped
cut -f 10 NAT.genes.all.bed | awk -F " " '{print $2}'| sort | uniq | sed 's~;~~g' > NAT.ids.txt

# Generate a bed file of the entire transcript of exons that were found to be NATs
grep -Ff NAT.ids.txt lincRNA.noOTs.bed > NAT.all.bed

#Create a list of all NAT transcripts
cut -f 10 NAT.all.bed | awk -F " " '{for(i=1;i<=NF;i++){if ($i ~/TCONS/) {print $i ".gene=" $2}}}'| sort | uniq | sed 's~;~~g' |sed 's~"~~g' | sed 's~zero_length_insertion=True~~g' |sed 's/^/>/' > NAT.all.txt

#create a bed file with all NATs and OTs removed
grep -vFf NAT.ids.txt lincRNA.noOTs.bed > lincRNA.postfilter.bed

# Move NATs to a new file
python /evolinc_docker/extract_sequences-1.py NAT.all.txt lincRNA.genes.fa Natural.antisense.transcripts.fa


ELAPSED_TIME_4=$(($SECONDS - $START_TIME_4))
echo "Elapsed time for Step 4 is" $ELAPSED_TIME_4 "seconds" >> ../$output/elapsed_time-evolinc-i.txt

# STEP 5: Generating final lincRNA.bed file and the final set of lincRNA sequences
START_TIME_5=$SECONDS
# Make a list from the lincRNA.postfilter.bed file in order to know which lincRNAs to extract from the lincRNA.genes.fa file
cut -f 10 lincRNA.postfilter.bed | awk -F " " '{for(i=1;i<=NF;i++){if ($i ~/TCONS/) {print $i ".gene=" $2}}}'| sort | uniq | sed 's~;~~g' |sed 's~"~~g' | sed 's~zero_length_insertion=True~~g' |sed 's/^/>/' > lincRNA.genes.filtered.uniq.genes

# Move lincRNA genes to a new file
python /evolinc_docker/extract_sequences-1.py lincRNA.genes.filtered.uniq.genes lincRNA.genes.fa All.lincRNAs.fa

# Sort the bed file for final output
sortBed -i lincRNA.postfilter.bed > lincRNA.bed

ELAPSED_TIME_5=$(($SECONDS - $START_TIME_5))
echo "Elapsed time for Step 5 is" $ELAPSED_TIME_5 "seconds" >> ../$output/elapsed_time-evolinc-i.txt

# STEP 6: Promoter extraction of filtered lincRNAs
START_TIME_6=$SECONDS
# Promoter extraction
sed 's/>//' lincRNA.genes.filtered.uniq.genes | cut -d "." -f 1 > lincRNA.genes.filtered.uniq.genes.mod

# Extract the coordinates from the cuffcompare file
grep -f lincRNA.genes.filtered.uniq.genes.mod ../$comparefile > lincRNA.genes.filtered.genes.gtf

# Extracting promoter coordinates from the gtf file for the transcripts
python /evolinc_docker/prepare_promoter_gtf.py lincRNA.genes.filtered.genes.gtf lincRNA.genes.filtered.genes.promoters.gtf

# Extracting fasta from the promoter coordinates
gffread lincRNA.genes.filtered.genes.promoters.gtf -w lincRNA.upstream.sequence.fa -g ../$referencegenome

ELAPSED_TIME_6=$(($SECONDS - $START_TIME_6))
echo "Elapsed time for Step 6 is" $ELAPSED_TIME_6 "seconds" >> ../$output/elapsed_time-evolinc-i.txt

# STEP 7: 
START_TIME_7=$SECONDS
# Update the cufflinks gtf file
python /evolinc_docker/update_gtf.py All.lincRNAs.fa ../$comparefile lincRNA.updated.gtf

ELAPSED_TIME_7=$(($SECONDS - $START_TIME_7))
echo "Elapsed time for Step 7 is" $ELAPSED_TIME_7 "seconds" >> ../$output/elapsed_time-evolinc-i.txt

# STEP 8:
START_TIME_8=$SECONDS
# Demographics for all lincRNA
quast.py All.lincRNAs.fa -R ../$referencegenome -G ../$referencegff -o lincRNA_demographics
sed 's/contig/lincRNA/g' lincRNA_demographics/report.txt > temp && mv temp lincRNA_demographics/report.txt
sed 's/contig/lincRNA/g' lincRNA_demographics/icarus.html > temp && mv temp lincRNA_demographics/icarus.html

# Identify the number of unique lincRNAs in bed file, and update the demographics file
uniquelincRNAcount=$(cut -f 10 lincRNA.bed | awk -F " " '{print $2}'| sort | uniq | grep -c "XLOC")
sed -i "4i Unique lincRNAs          $uniquelincRNAcount" lincRNA_demographics/report.txt
sed -i "s~# lincRNAs (>~# of total lincRNAs (including isoforms) (>~g" lincRNA_demographics/report.txt

# Demographics for Overlapping lincRNA
quast.py Overlapping.transcripts.fa -R ../$referencegenome -G ../$referencegff -o Overlapping.transcripts_demographics
sed 's/contig/lincRNA/g' Overlapping.transcripts_demographics/report.txt > temp && mv temp Overlapping.transcripts_demographics/report.txt
sed 's/contig/Overlapping.transcripts/g' Overlapping.transcripts_demographics/icarus.html > temp && mv temp Overlapping.transcripts_demographics/icarus.html

# Modify demographics file to include number of unique transcripts (unique XLOCs)
uniqueOTcount=$(cut -f 10 OT.all.bed | awk -F " " '{print $2}'| sort | uniq | grep -c "XLOC")
sed -i "4i Unique OT lncRNAs           $uniqueOTcount" Overlapping.transcripts_demographics/report.txt

# Change the terminology to reflect OTs
sed -i "s~# lincRNAs (>~# of total OT lncRNAs (including isoforms) (>~g" Overlapping.transcripts_demographics/report.txt
sed -i "s~lincRNA~OT lncRNA~g" Overlapping.transcripts_demographics/report.txt

# Demographics for natural antisense lincRNA
quast.py Natural.antisense.transcripts.fa -R ../$referencegenome -G ../$referencegff -o Natural.antisense.transcripts_demographics
sed 's/contig/lincRNA/g' Natural.antisense.transcripts_demographics/report.txt > Natural.antisense.transcripts_demographics/report.txt
sed 's/contig/Natural.antisense.transcripts/g' Natural.antisense.transcripts_demographics/icarus.html > Natural.antisense.transcripts_demographics/icarus.html

# Modify demographics file to include number of unique transcripts (unique XLOCs)
uniqueNATcount=$(cut -f 10 NAT.all.bed | awk -F " " '{print $2}'| sort | uniq | grep -c "XLOC")
sed -i "4i Unique NAT lncRNAs        $uniqueNATcount" Natural.antisense.transcripts_demographics/report.txt

# Change the terminology to reflect NATs
sed -i "s~# lincRNAs (>~# of total NAT lncRNAs (including isoforms) (>~g" Natural.antisense.transcripts_demographics/report.txt
sed -i "s~lincRNA~NAT lncRNA~g" Natural.antisense.transcripts_demographics/report.txt

ELAPSED_TIME_8=$(($SECONDS - $START_TIME_8))
echo "Elapsed time for Step 8 is" $ELAPSED_TIME_8 "seconds" >> ../$output/elapsed_time-evolinc-i.txt

# STEP 9:
START_TIME_9=$SECONDS
# Copy the files to the outputfiles
cp -r lincRNA.bed All.lincRNAs.fa lincRNA.upstream.sequence.fa lincRNA_demographics lincRNA.updated.gtf Overlapping.transcripts.fa Overlapping.transcripts_demographics Natural.antisense.transcripts.fa Natural.antisense.transcripts_demographics ../$output

ELAPSED_TIME_9=$(($SECONDS - $START_TIME_9))
echo "Elapsed time for Step 9 is" $ELAPSED_TIME_9 "seconds" >> ../$output/elapsed_time-evolinc-i.txt

# Optional STEP - 1:
START_TIME_O1=$SECONDS
# CAGE data
if [ ! -z $cagefile ]; then
     gff2bed < ../$cagefile > AnnotatedPEATPeaks.bed &&
     sortBed -i AnnotatedPEATPeaks.bed > AnnotatedPEATPeaks.sorted.bed &&
     closestBed -a lincRNA.bed -b AnnotatedPEATPeaks.sorted.bed -s -D a > closest_output.txt &&       
     python /evolinc_docker/closet_bed_compare.py closest_output.txt All.lincRNAs.fa lincRNAs.with.CAGE.support.annotated.fa &&
     cp lincRNAs.with.CAGE.support.annotated.fa ../$output
fi

ELAPSED_TIME_O1=$(($SECONDS - $START_TIME_O1))
echo "Elapsed time for Optional Step 1 is" $ELAPSED_TIME_O1 "seconds" >> ../$output/elapsed_time-evolinc-i.txt

# Optional STEP - 2:
START_TIME_O2=$SECONDS
# Known lincRNA
if [ ! -z $knownlinc ]; then
     gff2bed < ../$knownlinc > Atha_known_lncRNAs.bed &&
     sortBed -i Atha_known_lncRNAs.bed > Atha_known_lncRNAs.sorted.bed &&
     intersectBed -a lincRNA.bed -b Atha_known_lncRNAs.sorted.bed > intersect_output.txt &&
     python /evolinc_docker/interesect_bed_compare.py intersect_output.txt All.lincRNAs.fa lincRNAs.overlapping.known.lincs.fa &&
     cp lincRNAs.overlapping.known.lincs.fa ../$output
fi

ELAPSED_TIME_O2=$(($SECONDS - $START_TIME_O2))
echo "Elapsed time for Optional Step 2 is" $ELAPSED_TIME_O2 "seconds" >> ../$output/elapsed_time-evolinc-i.txt

# Pie chart if both CAGE and Knownlinc are given
if [ ! -z $cagefile ] && [ ! -z $knownlinc ] ; then
   python /evolinc_docker/lincRNA_fig.py All.lincRNAs.fa lincRNAs.with.CAGE.support.annotated.fa lincRNAs.overlapping.known.lincs.fa &&
   cp lincRNA_piechart.png ../$output
fi

# remove all the other files
rm -r ../transcripts_u_filter.fa.transdecoder_dir
rm ../*.fa*.*

echo "All necessary files written to" $output
echo "Finished Evolinc-part-I!"
echo `date`

ELAPSED_TIME=$(($SECONDS - $START_TIME))
echo "Total elapsed time is" $ELAPSED_TIME "seconds" >> ../$output/elapsed_time-evolinc-i.txt