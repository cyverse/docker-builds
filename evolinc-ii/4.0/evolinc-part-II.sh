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

while getopts ":hb:q:l:i:s:n:o:v:t:" opt; do
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
    o)
     output=$OPTARG # Output file
     ;;  
    t)
     species_tree=$OPTARG # Species tree
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
cat *reciprocal_sorted.gff $query_species.$query_species.coords.gff >All_Orthologs.gff
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

### Starting Phylogenetic steps ###
if [ ! -z $species_tree ];
then

  # Starting alignments
  echo "Starting alignments"
  ls * > alignment_list.txt
  sed -i '/alignment_list.txt/d' alignment_list.txt
  mkdir -p Final_results
  perl /Batch_MAFFT.pl alignment_list.txt
  echo "Finished with alignments, preparing files for RAxML if that option was selected"
  cd Final_results
  #If the number of fasta files in the current folder is equal to one proceed
  filenum=$(ls *.fasta | wc -l )
  if [ $filenum -lt 2 ]; then
  #The below code counts how many > are found in the single fasta file in the current directory
	  speciescount=$(grep -c ">" *.fasta)
	  #If the number of sequences in the file are less than 4, then do not create an aligned_list file
	  if [ $speciescount -lt 4 ]; then
		echo "not enough taxa represented to perform phylogenetics"
		else
		ls *.fasta >aligned_list.txt
		echo "Sufficient taxa represended for a single query lncRNA to perform phylogenetics"
	  fi
  else
  #If there are multiple files in the directory, the below script asks how many > each of the files has, only printing out the file names of those with more than 3 >s
	  grep -c ">" *.fasta | grep -v ':0$' | grep -v ':1$' | grep -v ':2$' | grep -v ':3$' | awk -F':' '{print $1}' >aligned_list.txt
	  #The below script is checking to make sure that the files have sufficient taxa
	  speciescount_2=$(grep -c ">" *.fasta | grep -v ':0$' | grep -v ':1$' | grep -v ':2$' | grep -v ':3$' | awk -F':' '{print $1}'|grep -c ".fasta")
	  if [ $speciescount_2 -gt 0 ]; then
	  echo "Sufficient taxa represented for multiple query lncRNAs to perform phylogenetics"
	  else
	  echo "Multiple files, however none with enough taxa present for phylogenetics"
	  fi
  fi

  #RAxML step

  echo "starting tree building"
  mkdir -p ../RAxML_families 
  perl /Batch_RAxML.pl aligned_list.txt
  mv aligned_list.txt ../RAxML_families
  cd ../RAxML_families
  mkdir Reconciled_trees
  mv ../Final_results/*.reduced Reconciled_trees
  perl /Batch_NOTUNG.pl aligned_list.txt ../../../$species_tree
  mv *png Reconciled_trees
  mv Reconciled_trees  ../../../$output
  cd ../Final_results

  # Generating summary of aligned lincRNA
  echo "Generating summary of aligned linRNA"
  python /Family_division_and_summary.py ../../../$sp_list
  grep -v "aligned_list.txt" summary.txt > final_summary.txt
  Rscript /final_summary.R -q ../../../$sp_list
  rm summary.txt
  mv lincRNA_barplot.png ../../../$output
  cd ../../
  Rscript /final_summary_table_gen.R -s ../$sp_list
  echo "Finished summary"
  
### End of phylogenetic step, the below code is if there is no phylogenetic step picked

else
  echo "Generating summary of aligned linRNA"
  python /Family_division_and_summary.py ../../$sp_list
  grep -v "aligned_list.txt" summary.txt > final_summary.txt
  cp final_summary.txt ../../$output/Instance_count.txt
  Rscript /final_summary.R -q ../../$sp_list
  cp -r */*.fasta .
  rm summary.txt
  mv lincRNA_barplot.png ../../$output
  cd ../
  Rscript /final_summary_table_gen.R -s ../$sp_list
  echo "Finished summary"
fi

# Modifying the final summary table to add ID's of knonw lincRNA and ID's of knonw sense and known antisense ID's

if compgen -G "../Homology_Search/*tested.out" > /dev/null && compgen -G "../Homology_Search/*mod.annotation.*sense.gff" > /dev/null;
then
        for i in ../Homology_Search/*tested.out; do
                if [[ -s $i ]] ; then
                   python /filter_lincRNA_sequences_annotation2.py "$i" "$i".mod
                   sed 's/_TBH_1//g' "$i".mod > temp && mv temp "$i".mod
                   python /filter_lincRNA_sequences_annotation3.py "$i".mod final_summary_table.tsv "$i".mod.sp.csv
                fi
        done
        for i in ../Homology_Search/*mod.annotation.*sense.gff; do
                if [[ -s $i ]] ; then
                   sed 's/_TBH_1//g' "$i" > temp && mv temp "$i"
                fi
        done
        Rscript /final_summary_table_all.R
        rm final_summary_table.tsv final_summary_table.mod.tsv
        mv final_summary_table.mod2.tsv final_summary_table.tsv
        mv final_summary_table.tsv ../$output

elif compgen -G "../Homology_Search/*tested.out" > /dev/null && ! compgen -G  "../Homology_Search/*mod.annotation.*sense.gff" > /dev/null
then
        for i in ../Homology_Search/*tested.out; do
                if [[ -s $i ]] ; then
                   python /filter_lincRNA_sequences_annotation2.py "$i" "$i".mod
                   sed 's/_TBH_1//g' "$i".mod > temp && mv temp "$i".mod
                   python /filter_lincRNA_sequences_annotation3.py "$i".mod final_summary_table.tsv "$i".mod.sp.csv
                fi
        done    
        Rscript /final_summary_table_all.R       
        rm final_summary_table.tsv
        mv final_summary_table.mod.tsv final_summary_table.tsv
        mv final_summary_table.tsv ../$output


elif ! compgen -G "../Homology_Search/*tested.out" > /dev/null && compgen -G  "../Homology_Search/*mod.annotation.*sense.gff" > /dev/null
then
   
        for i in ../Homology_Search/*mod.annotation.*sense.gff; do
                if [[ -s $i ]] ; then
                   sed 's/_TBH_1//g' "$i" > temp && mv temp "$i"
                fi
        done     
        Rscript /final_summary_table_all.R
        rm final_summary_table.tsv
        mv final_summary_table.mod.tsv final_summary_table.tsv
        mv final_summary_table.tsv ../$output

elif ! compgen -G "../Homology_Search/*tested.out" > /dev/null && ! compgen -G  "../Homology_Search/*mod.annotation.*sense.gff" > /dev/null
then
       
       mv final_summary_table.tsv ../$output  

fi


# Moving all the files to the output folder
cd ../
mv Homology_Search $output
mv BLAST_DB $output
mv Reciprocal_BLAST $output
mv Orthologs $output
rm stderr.out sample-instance-report.txt

echo "All necessary files written to" $output
echo "Finished Evolinc-part-II!"
echo `date`

ELAPSED_TIME=$(($SECONDS - $START_TIME))
echo "Total elapsed time is" $ELAPSED_TIME "seconds" >> $output/elapsed_time.txt
# END
