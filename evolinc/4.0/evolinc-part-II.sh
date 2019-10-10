#!/bin/bash
#author: Andrew Nelson; andrew.d.l.nelson@gmail.com
# Script to perform comparative genomic analysis of lincRNAs

usage() {
      echo ""
      echo "Usage : sh $0 -b BLASTing_list -l Query_lincRNA -q Query_species -i Input_folder -s Species_list -v evalue -t species_tree -o Output_folder"
      echo ""

cat <<'EOF'
  -b </path/to/BLASTing_list in tab-delimited format>

  -l </path/to/query lincRNA>

  -q <query species>

  -i </path/to/input folder>

  -s </path/to/species list>

  -v <e-value>
  
  -t </path/to/species/tree>

  -o </path/to/output folder>

  -h Show this usage information

EOF
    exit 0
}

while getopts ":hb:q:l:i:s:o:v:t" opt; do
  case $opt in
    b)
     Blasting_list=$OPTARG #This is a five-column tab-delimited list in the following order:
		# subject_genome lincRNA_fasta	query_species (four letter Genus-species designation ie., Gspe)	subject_species (same four letter abbreviation as subject_species)	               subject_gff (in fasta_format) # All of these files should be in the current working folder
     ;;
    l)
    query_lincRNA=$OPTARG # Query lincRNA file here
     ;;
    q)
    query_species=$OPTARG # This is entirely left to the user
     ;;
    i)
    input_folder=$OPTARG # Input folder
     ;;
    s)
     sp_list=$OPTARG # Species list 
     ;;
    v)
     value=$OPTARG # e-value
     ;; 
    t)
     species_tree=$OPTARG
     ;;
    o)
     output=$OPTARG # Output file
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

echo "BEGIN!"
echo `date`

START_TIME=$SECONDS

# Make all necessary folders
mkdir $output

mkdir BLAST_DB
mkdir Homology_Search
mkdir Reciprocal_BLAST
mkdir Orthologs

# Make lincRNA list from the user provide lincRNA fasta file
grep ">" $query_lincRNA | sed 's/>//' > Orthologs/lincRNA.list

# Initiate search for putative orthologs
echo "***Starting lincRNA to Genome Comparisons***"
python /startup_script.py -b $Blasting_list -i $input_folder -v $value
echo "***Finished with lincRNA to Genome Comparisons***"

#Create a list of all genomes, lincRNA ortholog, and gff files to set up reciprocal BLAST
cd Reciprocal_BLAST
find . -maxdepth 1 -name "*genome*" >Reciprocal_chrom_list.txt
find . -maxdepth 1 -name "*orthologs*" >Reciprocal_lincRNAs_list.txt
find . -maxdepth 1 -name "*coords*" >Reciprocal_lincRNAs_coord_list.txt
sed -i 's~./~~g' Reciprocal_chrom_list.txt
sed -i 's~./~~g' Reciprocal_lincRNAs_list.txt
sed -i 's~./~~g' Reciprocal_lincRNAs_coord_list.txt
sort -n Reciprocal_chrom_list.txt -o Reciprocal_chrom_list.txt
sort -n Reciprocal_lincRNAs_list.txt -o Reciprocal_lincRNAs_list.txt
sort -n Reciprocal_lincRNAs_coord_list.txt -o Reciprocal_lincRNAs_coord_list.txt
:|paste Reciprocal_chrom_list.txt - Reciprocal_lincRNAs_list.txt - Reciprocal_lincRNAs_coord_list.txt >Reciprocal_list.txt
sed -i 's~\t\t~\t~g' Reciprocal_list.txt
sed -i "s/$/\t$query_species.$query_species.coords.gff/" Reciprocal_list.txt
sed -i "s/$/\t$query_species/" Reciprocal_list.txt
sed -i "s/$/\t$query_species.genome.fasta/" Reciprocal_list.txt
#sed -i "s~\.putative~_putative~g" Reciprocal_list.txt
sed -i "s~\.genome~_genome~g" Reciprocal_list.txt
rm Reciprocal_chrom_list.txt
rm Reciprocal_lincRNAs_list.txt
rm Reciprocal_lincRNAs_coord_list.txt

# Confirm the reciprocity of the putative orthologs
echo "***Starting Reciprocal Search***"
python /Reciprocal_BLAST_startup_script.py -b Reciprocal_list.txt -v $value
echo "***Finished with Reciprocal Search***"

#Create a CoGe viewable bed file of TBH
cat *TBH.out.gff $query_species.$query_species.coords.gff >All_Orthologs.gff
awk -F '\t' '{print $1 "\t" $4 "\t" $5 "\t" $6 "\t" $8 "\t" $7 "\t" $2 "\t" $3 "\t" $6 "\t" $9}' All_Orthologs.gff >All_orthologs.bed
sort -k 1,1 -k 2,2n All_orthologs.bed > ../$output/All_orthologs_for_viewing.bed

#Starting building families
echo "***Creating Families of similar sequences***"
cd ../Orthologs
cat *.fasta > All_orthologs.fasta
perl /find_from_list.pl lincRNA.list All_orthologs.fasta
cd lincRNA_families
for i in *FASTA; do mv $i "`basename $i _.FASTA`.fasta"; done
echo "Finished Creating Families of similar sequences"

# Starting alignments
echo "Starting alignments"
ls * > alignment_list.txt
sed -i '/alignment_list.txt/d' alignment_list.txt
mkdir -p Final_results
perl /Batch_MAFFT.pl alignment_list.txt
cd Final_results
grep -c ">" * | grep -v ':0$' | grep -v ':1$' | grep -v ':2$' | grep -v ':3$' | awk -F':' '{print $1}' > aligned_list.txt
echo "Finished alignments"
echo "Finished alignments, starting tree building"
#ls * >aligned_list.txt
#sed -i '/aligned_list.txt/d' aligned_list.txt
# The species_tree isn't actually required here, but I am using it as a placeholder to tell the program if the user wants to proceed to the RAxML step.
# Once we have a different activation signal, we can get rid of this file requirement.

#RAxML optional step
if [ ! -z $species_tree ]; 
then 
  mkdir -p ../RAxML_families
  perl /Batch_RAxML.pl aligned_list.txt
  rm aligned_list.txt
else
  rm aligned_list.txt
fi

# Generating summary of aligned lincRNA
echo "Generating summary of aligned linRNA"
python /Family_division_and_summary.py ../../../$sp_list
grep -v "aligned_list.txt" summary.txt > final_summary.txt
Rscript /final_summary.R -q ../../../$sp_list
rm summary.txt
mv lincRNA_barplot.png ../../../
cd ../../
Rscript /final_summary_table_gen.R -s ../$sp_list
#mv lincRNA_families/Final_results ../../$output
mv final_summary_table.tsv ../
echo "Finished summary"

# Moving all the files to the output folder
cd ../
mv Homology_Search $output
mv BLAST_DB $output
mv Reciprocal_BLAST $output
mv Orthologs $output
mv stderr.out sample-instance-report.txt $output

echo "All necessary files written to" $output
echo "Finished Evolinc-part-II!"
echo `date`

ELAPSED_TIME=$(($SECONDS - $START_TIME))
echo "Total elapsed time is" $ELAPSED_TIME "seconds" >> $output/elapsed_time.txt
# END
