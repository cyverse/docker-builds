#!/bin/bash 
set -x

#This script first sets up ontologizer and downloads the latest .obo and assoication files
#The user can enter in a taxon id and the script will download the corresponding gene association file from www.geneontology.org 
#If the gene association file cannot be found in geneontology.org this tool creates a GAF for the taxon using QuickGO.
#For this script to function properly subversion (svn) needs to be installed!
#Author: Kent Guerriero
#Modified January 2020 to work with new Gene Ontology links and file names.
#Global vars
taxonid=0 #This is the taxonID from the user
g_arg="not selected" #If the user has a GAF they can use their own
p_arg="not selected" #The user needs to enter this (this is their population set)
o_arg="not selected" #This is the location where the output file is written
geneOboFileName="go-basic.obo"

#This is the name of the filtered files on geneOntology.org that this script will try to use first
#I don't know of a practical way of getting these file names without hardcoding it...
geneOntologyFileNames="Arabidopsis thaliana     3702   :tair.gaf.gz,
Aspergillus nidulans 162425             :aspgd.gaf.gz,
Bos taurus      9913    :goa_cow.gaf.gz,
Caenorhabditis elegans  6239    :wb.gaf.gz,
Candida albicans        5476    :cgd.gaf.gz,
Canis lupus familiaris  9616    :goa_dog.gaf.gz,
Danio rerio     7955    :zfin.gaf.gz,
Dictyostelium discoideum        44689   :dictybase.gaf.gz,
Drosophila melanogaster 7227    :fb.gaf.gz,
Escherichia coli        562     :ecocyc.gaf.gz,
Gallus gallus   9031    :goa_chicken.gaf.gz,
Homo sapiens    9606    :goa_human.gaf.gz,
Leishmania major        5664    :genedb_lmajor.gaf.gz,
Mus musculus    10090   :mgi.gaf.gz,
Pseudomonas aeruginosa 3716	:pseudocap.gaf.gz,
Rattus norvegicus       10116   :rgd.gaf.gz,
Saccharomyces cerevisiae        4932    :sgd.gaf.gz,
Schizosaccharomyces pombe       4896    :pombase.gaf.gz,
Sus scrofa      9823    :goa_pig.gaf.gz,
Trypanosoma brucei      5691    :genedb_tbrucei.gaf.gz,"


#checks that the option to select a taxon is set
if [ "$1" == "-t" ];then
   taxonid=$2
   echo "User entered $taxonid"
elif [ "$1" == "--userGAF" ];then #if not using taxon id, user has a gaf they want to use 
   taxonid=-1
   fileName=$2
   echo "User is using their own GAF file located at $2"
fi

#checks to see that the user entered a taxonid or a GAF file
if [ $taxonid == '0' ];then
   exit 1
fi

#This shifts the current number of arguments by 2 (removes arguments specific to this script)
#The remaining arguments are passed into ontologizer.jar
shift 2

if [ "$1" == "--userGAF" ];then #if user has entered both --userGAF and -t, just use userGAF 
   taxonid=-1
   fileName=$2
   echo "User is using their own GAF file located at $2"
   shift 2
fi
ontologizerArgs="$@"


#Check to see if .obo file is up to date(using timestamp)
#If not download the latest .obo file
   downloadLoc="http://current.geneontology.org/ontology/"
   downloadLoc=$downloadLoc""$geneOboFileName
   wget $downloadLoc

#Checks if .obo downloaded successfully
if [ $? -ne 0 ];then 
   echo "There was a problem downloading the .obo file"
   exit 1
fi


#Generic help statment for the script
#This explains what ontologizer does/ needs and what this script does/needs
printHelpStatement(){
   echo "
This tool will download the latest gene association files (.gaf) given a taxon id. 
If the GAF is a multispecies file, this program will filter out the non relevant species.
This tool also updates the GeneOntology file(.obo) and the Ontologizer program (if needed). 


Ontologizer requires 4 input files to run correctly 
(1) GenoOntology (.obo) file
   We download the latest version for you
(2) The Gene Association file
   This is what we search for when you enter your taxon ID
   If we cannot find this, you will need to provide this to Ontologizer.
(3) The population set 
   You will need to provide this to Ontologizer.
(4) The study set
   You will need to provide this to Ontologizer
" 
}


#This function gets the taxonid from the user 
#It then searches the geneOFileLocations.txt file to see if the file is located at geneOntology.org
searchWithTaxonId(){
   taxIdLocation=""
   #First we will use user input to look for taxon ID in files at geneOntology.org
   taxIdLocation=$(echo $geneOntologyFileNames | sed "s/,/,\n/g" | egrep " $taxonid ") 
   searchForFile
}

#Creates the GAF file from QuickGo
createFileFromQuickGOSource(){
   echo "File is being created, please wait"
   fileName="gene_association$taxonid.gaf" 
   wget -O ./$fileName "http://www.ebi.ac.uk/QuickGO/GAnnotation?format=gaf&limit=-1&tax=$taxonid"
}

#Creates the GAF from GeneOntology.org
createFileFromGeneOSource(){
   downloadLocation="http://current.geneontology.org/annotations/"
   downloadLocation=$downloadLocation""$fileName
   wget $downloadLocation
   echo "Found file, downloading..."
   
#Checks to see if there was a problem downloading the gene association file 
   if [ $? -ne 0 ]; then
      echo "There was a problem downloading the file, you may need to manually download it at geneontology.org"
      exit 1
   else
#Unzips and moves the GAF to the directory where the script is
      gunzip $fileName
      echo "File saved in $( echo $PWD ) "
   fi
}

#This function ensures that the file found from the user given taxon ID is correct. 
#The function then downloads the correct GAF files according to the taxon ID 
#If we cannot locate the taxon id in the files from geneOntology.org we will use QuickGO to create the GAF
searchForFile()
{
   isInfoCorrect='n'
#Checks to see if the file found is the correct file (asks user)
   if [ "$taxIdLocation" != "" ]; then
      echo "FOUND: $taxIdLocation"
      isInfoCorrect='y'
      fileName=$(echo $taxIdLocation | egrep -o ":.*," | tr -d ':'| tr -d ',')
      fileLocation=$fileName
      createFileFromGeneOSource
#This part of the code checks for a second possible GAF
      fileName=$(echo $taxIdLocation | egrep -o ",.*" | tr -d ':'| tr -d ',')  
      if [ "$fileName" != "" ]; then
         createFileFromGeneOSource
      else
         fileName=$fileLocation
      fi
      fileName=$(echo $fileName | sed 's/.gz//g')
   else  #If we cant find the taxonid in the files at geneontology.org we use QuickGO to generate the GAF
      echo "Taxon id $taxonid NOT FOUND in GAF file from GO.org. Creating GAF file using quickGO"
         createFileFromQuickGOSource     
   fi
}

#Display welcome message
#This will only display if run via the command line
echo "
Welcome to Ontologizer on iPlant
This script is designed to download and update files assoicated with Ontologizer
"
#TaxonID will be -1 if the user gave their own GAF file
if [ $taxonid != "-1" ]; then
   searchWithTaxonId
fi

#This calls Ontologizer with the proper arguments
#Note that this is modified to work on the iPlant DE taking the taxon id or the location of the user GAF first (-t or --userGAF). If this order is modified this script may not function correctly.
/usr/bin/java -jar /ontologizer/Ontologizer.jar $ontologizerArgs -a $fileName -g ./$geneOboFileName

rm $fileName $geneOboFileName

exit 0
