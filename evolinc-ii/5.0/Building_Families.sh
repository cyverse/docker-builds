#!/bin/bash
#author: Andrew Nelson; andrew.d.l.nelson@gmail.com
# Script to process cuffcompare output file to generate lincRNA
# Usage: 
# sh Building_Families.sh  -g subjectgenome.fa -s subject_species -q query_species -l lincRNAs.fa -e subject_gff -k known_lincRNAs -v e_value

while getopts ":l:q:s:f:k:e:g:" opt; do
  case $opt in
    l)
      lincRNAfasta=$OPTARG
    ;;
    q)
      query_species=$OPTARG
      ;;
    s)
      subject_species=$OPTARG
      ;;
    f)
      subject_gff=$OPTARG
      ;;
    g)
      subject_genome=$OPTARG
      ;;
    k)
      known_lincRNAs=$OPTARG
      ;;
    e)
      value=$OPTARG
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

echo "Building_Families.sh script is reporting this"
echo "Starting! on $subject_species"

# Formatting the genome
makeblastdb -logfile stderr.out -in $subject_genome -dbtype nucl -out BLAST_DB/$subject_genome.blast.out
### Check to see if query FASTA file contains only one lincRNA
count=$(grep -c ">" $lincRNAfasta)
param=2
if [ "$count" -ge "$param" ]; then
   touch $lincRNAfasta
else
   cat $lincRNAfasta $lincRNAfasta >temp && mv temp $lincRNAfasta
   sed -i 's~.>~\n>~g' $lincRNAfasta
fi

# Blasting the lincRNA transcripts against the genome to find out the location on the genome And identify if there are any paralogs in the query genome.
blastn -logfile stderr.out -query $lincRNAfasta -db BLAST_DB/$subject_genome.blast.out -num_threads 4 -penalty -2 -reward 1 -gapopen 5 -gapextend 2 -dust no -word_size 8 -evalue $value -outfmt "6 qseqid sseqid pident length qlen qstart qend sstart send evalue bitscore" -out Homology_Search/$subject_species.out

# Remove spaces in the blastout files
sed 's/ //g' Homology_Search/$subject_species.out > Homology_Search/$subject_species.stripped.out

# Convert blast result to gff
perl /blast2gff.pl -i Homology_Search/$subject_species.stripped.out -s $subject_species -o Homology_Search/$subject_species.out.gff

###Sort the gff file and merge the start and stop including the intermediate sequences. Check the merge_close_hits.py script. No need to do this since the script is doing exactly what it is supposed to do.. We can alter the BLAST2GFF.PL file to order the columns in the way that gffread likes. Look at the code starting at line #175 for reordering.
# Merge close hits in the gff file # Merge the start and stop coordinates
python /merge_close_hits.py Homology_Search/$subject_species.out.gff Homology_Search/$subject_species.out.merged.gff

grep "TBH" Homology_Search/$subject_species.out.merged.gff >Homology_Search/$subject_species.out.TBH.only.gff
#The below lines grab five of the top non-TBH hits (in case of large, multi-hit families we don't want to grab too many). This is informative for examining duplication events.
grep -v "TBH" Homology_Search/$subject_species.out.merged.gff >Homology_Search/$subject_species.out.no_TBH_precursor.only.gff
cut -f 1 Homology_Search/$subject_species.out.no_TBH_precursor.only.gff | grep '._1$\|._2$\|._3$\|._4$\|._5$' | sed 's~$~\t~g' >Homology_Search/$subject_species.out.no_TBH.only.gff.list.txt
grep -f Homology_Search/$subject_species.out.no_TBH.only.gff.list.txt Homology_Search/$subject_species.out.no_TBH_precursor.only.gff >Homology_Search/$subject_species.out.no_TBH.only.gff

# Convert gtf to fasta
# Change the file format since the gffread requires the chromosome to be first column
awk '{print $2 "\t" $1 "\t" $3 "\t" $4 "\t" $5 "\t" $6 "\t" $7 "\t" $8 "\t" $9 " " $10  " " $11 " " $12}' Homology_Search/$subject_species.out.merged.gff > temp && mv temp Homology_Search/$subject_species.out.merged.gff
awk '{print $2 "\t" $1 "\t" $3 "\t" $4 "\t" $5 "\t" $6 "\t" $7 "\t" $8 "\t" $9 " " $10  " " $11 " " $12}' Homology_Search/$subject_species.out.TBH.only.gff > temp && mv temp Homology_Search/$subject_species.out.TBH.only.gff
awk '{print $2 "\t" $1 "\t" $3 "\t" $4 "\t" $5 "\t" $6 "\t" $7 "\t" $8 "\t" $9 " " $10  " " $11 " " $12}' Homology_Search/$subject_species.out.no_TBH.only.gff > temp && mv temp Homology_Search/$subject_species.out.no_TBH.only.gff
gffread Homology_Search/$subject_species.out.merged.gff -g $subject_genome -w Homology_Search/$subject_species.$query_species.orthologs.fasta
gffread Homology_Search/$subject_species.out.TBH.only.gff -g $subject_genome -w Homology_Search/$subject_species.$query_species.TBH.orthologs.fasta
gffread Homology_Search/$subject_species.out.no_TBH.only.gff -g $subject_genome -w Homology_Search/$subject_species.$query_species.no_TBH.orthologs.fasta

# Optional argument - 1 (transcriptome) 
# Here we can use bedtools to compare the Hsap.out.merged.gff (lincRNA) against the $subject_gff. the code will run something like so: 
# Bedtools can check for overlap in either strand. To work with downstream scripts, we can run bedtools twice, one checking for sense and making a list, and the other checking for antisense #(and making a list).
if [ ! -z $subject_gff ]; then 
    # sense
    intersectBed -a Homology_Search/$subject_species.out.TBH.only.gff -b $subject_gff -f 0.1 -u -s > Homology_Search/$subject_species.annotation.sense.gff &&
    #the -f of 0.1 indicates a 10% overlap. This may be too much or not enough. We will have to check emperically. -u reports the mere presence of any unique intersecting gene.
    # antisense
    # This is for capturing the GeneID of the overlapping gene
    intersectBed -wb -a Homology_Search/$subject_species.out.TBH.only.gff -b $subject_gff -f 0.1 -s > Homology_Search/$subject_species.mod.annotation.sense.gff &&
    awk '{print $2}' Homology_Search/$subject_species.annotation.sense.gff > Homology_Search/$subject_species.annotation.sense.list.txt &&
    # Retain only the second column (gene IDs) and save as a list
    intersectBed -a Homology_Search/$subject_species.out.TBH.only.gff -b $subject_gff -f 0.1 -u -S > Homology_Search/$subject_species.annotation.antisense.gff &&
    # This is the same as above, except for -S searches only for those overlapping candidates that are on the antisense strand.
    # This is for capturing the GeneID of the overlapping gene
    intersectBed -wb -a Homology_Search/$subject_species.out.TBH.only.gff -b $subject_gff -f 0.1 -S > Homology_Search/$subject_species.mod.annotation.antisense.gff &&
    awk '{print $2}' Homology_Search/$subject_species.annotation.antisense.gff > Homology_Search/$subject_species.annotation.antisense.list.txt &&
    # Retain only the second column (gene IDs) and save as a list
    # assign annotation
    # sense
    python /assign_sense_annotation.py Homology_Search/$subject_species.$query_species.TBH.orthologs.fasta Homology_Search/$subject_species.annotation.sense.list.txt Homology_Search/$subject_species.$query_species.orthologs.sense.renamed.fasta &&
    # assign_sense_annotation.py works as .py file_with_sequences_to_rename.fasta list_of_headers_to_rename >renamed_output.fasta
    #sense + antisense
    python /assign_antisense_annotation.py Homology_Search/$subject_species.$query_species.orthologs.sense.renamed.fasta Homology_Search/$subject_species.annotation.antisense.list.txt Homology_Search/$subject_species.$query_species.orthologs.renamed.fasta
else
		mv Homology_Search/$subject_species.$query_species.TBH.orthologs.fasta Homology_Search/$subject_species.$query_species.orthologs.renamed.fasta
fi

# Optional argument - 2 (Known lincRNAs)
if [ ! -z $known_lincRNAs ]; then 
    # Formatting the known lincRNA's
    makeblastdb -logfile stderr.out -in $known_lincRNAs -dbtype nucl -out BLAST_DB/$known_lincRNAs.blast.out &&
    # Blasting to known lincRNAs 
    blastn -logfile stderr.out -query Homology_Search/$subject_species.$query_species.orthologs.renamed.fasta -db BLAST_DB/$known_lincRNAs.blast.out -num_threads 2 -penalty -2 -reward 1 -gapopen 5 -gapextend 2 -dust no -word_size 8 -evalue $value -outfmt "6 qseqid sseqid pident length qlen qstart qend sstart send evalue bitscore" -out Homology_Search/$subject_species.$query_species.orthologs.renamed.lincRNAs_tested.out &&
    # Filtering the output
    python /filter_lincRNA_sequences_annotation.py Homology_Search/$subject_species.$query_species.orthologs.renamed.lincRNAs_tested.out Homology_Search/$subject_species.lincRNA_annotation.list.txt &&
    # Assign the annotation of lincRNA to the known lincRNA
    python /assign_annotation_lincRNA.py Homology_Search/$subject_species.$query_species.orthologs.renamed.fasta Homology_Search/$subject_species.lincRNA_annotation.list.txt Homology_Search/$subject_species.$query_species.orthologs.lincRNA_tested.renamed.fasta
else
	   mv Homology_Search/$subject_species.$query_species.orthologs.renamed.fasta Homology_Search/$subject_species.$query_species.orthologs.lincRNA_tested.renamed.fasta
fi

# Sort the fasta headers
#perl Sort_FASTA_Alp.pl -r Homology_Search/$subject_species.$query_species.orthologs.lincRNA_tested.renamed.fasta >Homology_Search/$subject_species.$query_species.orthologs_alpha.fasta

# Remove duplicate sequences
perl /Remove_dup_seqs.pl Homology_Search/$subject_species.$query_species.orthologs.lincRNA_tested.renamed.fasta

### Initiate reiterative BLAST here:
grep ">" Homology_Search/$subject_species.$query_species.orthologs.lincRNA_tested.renamed.fasta.dup_removed.fasta |sed 's~^......~~g'| sed 's~_Known_lincRNA~~g'| sed 's~_Known_Gene_Sense~~g' | sed 's~_Known_Gene_Antisense~~g' | sed 's~_TBH_1~~g' | sed 's~_.$~~g' | sed 's~_..$~~g' | sed 's~_...$~~g' | sed 's~>~~g' | sort -u > Homology_Search/List_of_identified_putative_orthologs.txt
grep ">" $lincRNAfasta > Homology_Search/List_of_all_query_lincRNAs.txt
# Search for lincRNA IDs that have not been found with the previous BLAST search. If non are found, still populate the file with the title "List_of_non_identified_lincRNAs" so the next step doesn't fail.
grep -vwf Homology_Search/List_of_identified_putative_orthologs.txt Homology_Search/List_of_all_query_lincRNAs.txt > Homology_Search/List_of_non_identified_query_lincRNAs.txt
sed -i 's~>~~g' Homology_Search/List_of_non_identified_query_lincRNAs.txt
perl /singleline.pl $lincRNAfasta > Homology_Search/query_lincRNAs_singleline.fasta
grep -A 1 -f Homology_Search/List_of_non_identified_query_lincRNAs.txt Homology_Search/query_lincRNAs_singleline.fasta | sed 's~--~~g' > Homology_Search/$query_species.$subject_species.non_identified_from_first_round.fasta
perl /split_fasta_into_200nt.pl Homology_Search/$query_species.$subject_species.non_identified_from_first_round.fasta > Homology_Search/$query_species.$subject_species.into_200_bin.fasta

##
##
##
#At this point, we are merely repeating the entire BLAST step with smaller pieces

# Blasting the lincRNA transcripts against the genome to find out the location on the genome And identify if there are any paralogs in the query genome.
blastn -logfile stderr.out -query Homology_Search/$query_species.$subject_species.into_200_bin.fasta -db BLAST_DB/$subject_genome.blast.out -num_threads 4 -evalue $value -outfmt "6 qseqid sseqid pident length qlen qstart qend sstart send evalue bitscore" -out Homology_Search/$subject_species.into_200_bin.fasta.out

# Remove spaces in the blastout files
sed 's/ //g' Homology_Search/$subject_species.into_200_bin.fasta.out > Homology_Search/$subject_species.into_200_bin_stripped.out

# Convert blast result to gff
perl /blast2gff.pl -i Homology_Search/$subject_species.into_200_bin_stripped.out -s $subject_species -o Homology_Search/$subject_species.into_200_bin.gff
skip=Homology_Search/$subject_species.into_200_bin.gff
if [[ $(du -h "$skip" | cut -f 1) = 0 ]]; then
  echo "No Matches Found in reiterative BLAST"
  echo ">No_further_matches_found_by_reiterative_BLAST" >Homology_Search/$subject_species.$query_species.orthologs.into_200_bin_lincRNA_tested.renamed.fasta.dup_removed.fasta
  echo "#No_further_matches_found_by_reiterative_BLAST" >Homology_Search/$subject_species.into_200_bin.TBH.only.gff
else
  echo "Reiterative BLAST found at least one hit, proceeding to rename and merge, if applicable"
fi
###Sort the gff file and merge the start and stop including the intermediate sequences. Check the merge_close_hits.py script. No need to do this since the script is doing exactly what it is supposed to do.. We can alter the BLAST2GFF.PL file to order the columns in the way that gffread likes. Look at the code starting at line #175 for reordering.
# Merge close hits in the gff file # Merge the start and stop coordinates
python /merge_close_hits.py Homology_Search/$subject_species.into_200_bin.gff Homology_Search/$subject_species.into_200_bin.merged.gff

grep "TBH" Homology_Search/$subject_species.into_200_bin.merged.gff > Homology_Search/$subject_species.into_200_bin.TBH.only.gff
# Convert gtf to fasta (Instead of using gtftocdna, we will use gffread)
# Change the file format since the gffread requires the chromosome to be first column
awk '{print $2 "\t" $1 "\t" $3 "\t" $4 "\t" $5 "\t" $6 "\t" $7 "\t" $8 "\t" $9 " " $10  " " $11 " " $12}' Homology_Search/$subject_species.out.merged.gff > temp && mv temp Homology_Search/$subject_species.out.merged.gff
awk '{print $2 "\t" $1 "\t" $3 "\t" $4 "\t" $5 "\t" $6 "\t" $7 "\t" $8 "\t" $9 " " $10  " " $11 " " $12}' Homology_Search/$subject_species.into_200_bin.TBH.only.gff > temp && mv temp Homology_Search/$subject_species.into_200_bin.TBH.only.gff
#gffread Homology_Search/$subject_species.out.merged.gff -g $subject_genome -w Homology_Search/$subject_species.$query_species.orthologs.fasta
gffread Homology_Search/$subject_species.into_200_bin.TBH.only.gff -g $subject_genome -w Homology_Search/$subject_species.$query_species.into_200_bin.TBH.orthologs.fasta
echo "Finished with reiterative BLAST"
cat Homology_Search/$subject_species.$query_species.into_200_bin.TBH.orthologs.fasta Homology_Search/$subject_species.$query_species.orthologs.lincRNA_tested.renamed.fasta.dup_removed.fasta > Homology_Search/$subject_species.$query_species.200nt_plus_full_length_searches.fasta
cat Homology_Search/$subject_species.into_200_bin.TBH.only.gff Homology_Search/$subject_species.out.TBH.only.gff > Homology_Search/$subject_species.200nt_plus_full_length_searches.gff
sort -k 1,1 -k 4,4n Homology_Search/$subject_species.200nt_plus_full_length_searches.gff > Homology_Search/$subject_species.200nt_plus_full_length_searches_sorted.gff
### End of reiterative BLAST section

cat Homology_Search/$subject_species.$query_species.200nt_plus_full_length_searches.fasta Homology_Search/$subject_species.$query_species.no_TBH.orthologs.fasta > Homology_Search/$subject_species.$query_species.200nt_plus_all_put_homologs.fasta

# Move genome and putative ortholog files for reciprocal BLAST to reciprocal BLAST folder. Also, rename putative orthologs file
cp $subject_genome Reciprocal_BLAST/
cp Homology_Search/$subject_species.$query_species.200nt_plus_full_length_searches.fasta Reciprocal_BLAST/$subject_species.$query_species.putative_orthologs.fasta
cp Homology_Search/$subject_species.200nt_plus_full_length_searches_sorted.gff Reciprocal_BLAST/$subject_species.$query_species.coords.gff
cp Homology_Search/$subject_species.$query_species.200nt_plus_all_put_homologs.fasta Orthologs/$subject_species.$query_species.putative_orthologs.fasta
echo "Finished with" $subject_species
