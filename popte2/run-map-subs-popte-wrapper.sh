#!/bin/bash

set -x
set -e

usage() {
      echo ""
      echo "Usage : sh $0 -g <TE-combined-reference> -i <Te-hierachy> -l <fastq1-reads> -r <fastq2-reads> -o <output_folder> -t num_threads -m <mincount> -c <targetcov> -q <mapqual> -d <mode> -x <maxdi> -n <mind> -a <maxd>"
      echo ""

cat <<'EOF'
  
  ###### Command line options ##########
  -g The repeat-masked reference genome
  -i Te-hierachy
  -l fastq file with first read pairs (gziped allowed)
  -r fastq file with second read pairs (gziped allowed)
  -o Output folder
  -t Number of threads
  -5 5' trim 
  -m the minimum count of a TE insertion; default=2.0
  -c the target coverage of the output file [int]	 
  -q minimum mapping quality; default=15
  -d joint|separate
  -x The maximum disagreement for the strand of the TE insertion fraction of reads
  -n The minimum distance between signatures; default=-100
  -a The maximum distance between signatures; default=500

EOF
    exit 0
}

while getopts ":hg:i:l:r:o:t:m:c:q:d:x:n:a:" opt; do
  case $opt in
    g)
    refg=$OPTARG # Reference genome file
     ;;
    i)
	hier=$OPTARG
	 ;;
    l)
    left_reads+=("$OPTARG") # Left reads
     ;;
    r)
    right_reads=("$OPTARG") # Right reads
     ;;
    o)
    of=$OPTARG # Samoutput file
     ;;
    t)
    threads=$OPTARG # Number of threads
     ;;
    m)
    mincount=$OPTARG 
     ;;
    c)
    targetcov=$OPTARG
     ;;
    q)
	mapqual=$OPTARG
	 ;;
	d)
	mode=$OPTARG
	 ;;
	x)
	maxdi=$OPTARG
	 ;;
	n)
	mind=$OPTARG
	 ;;
	a)
	maxd=$OPTARG
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


########################    Create directory     ############################################
mkdir $of

########################    map     #########################################################
bwa bwasw -t $threads $refg $left_reads |samtools view -Sb - > $of/map-read1.bam
bwa bwasw -t $threads $refg $right_reads |samtools view -Sb - > $of/map-read2.bam
java -jar popte2-v1.10.04.jar se2pe --fastq1 $left_reads --fastq2 $right_reads --bam1  $of/map-read1.bam --bam2 $of/map-read2.bam --sort --output $of/map-pe.sort.bam

#####################    POPTE2    ###########################################################
# ppileup
java -jar popte2-v1.10.04.jar ppileup --bam $of/map-pe.sort.bam --map-qual $mapqual --hier $hier --output $of/pp.gz

# subsample
java -jar popte2-v1.10.04.jar subsamplePpileup --ppileup $of/pp.gz --output $of/pp.subs.gz --target-coverage $targetcov

# rest
java -jar popte2-v1.10.04.jar identifySignatures --ppileup $of/pp.subs.gz --mode $mode --min-count $mincount --output $of/te.signatures
java -jar popte2-v1.10.04.jar updatestrand --bam $of/map-pe.sort.bam --signature $of/te.signatures --output $of/testrand.signatures --hier $hier --map-qual $mapqual --max-disagreement $maxdi
java -jar popte2-v1.10.04.jar frequency --ppileup $of/pp.subs.gz --signature $of/testrand.signatures --output $of/te.freqsignatures
java -jar popte2-v1.10.04.jar pairupsignatures --signature $of/te.freqsignatures --ref-genome $refg --hier $hier --min-distance $mind --max-distance $maxd --output $of/tes.finalresult.txt
