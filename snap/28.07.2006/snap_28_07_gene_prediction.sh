#!/bin/bash
# Upendra Kumar Devisetty
# Wrapper script snap tool
# 03-19-2016
# Usage: 
# sh evolinc-part-I.sh -c sample.data/cuffcompare_out_annot_no_annot.combined.gtf -g sample.data/Brassica_rapa_v1.2_genome.fa -r sample.data/Brassica_rapa_v1.2.gff -b sample.data/TE_RNA_transcripts.fa -t CAGE_file.txt -x 2> errors.txt > output.txt

usage() {
      echo ""
      echo "Usage : sh $0 -m <HMM file> -f <FASTA file> -g <gff file> -a <ace file> -p <plus strand> -e <negative strand> -n <name> -i <aminoacid file> -t <transcripts file> -l <lower case>"
      echo ""

cat <<'EOF'

  -m </path/to/HMM filee>

  -f </path/to/FASTA file>

  -g </path/to/gff file>

  -a </path/to/ace file>

  -p <predict on plus strand only>

  -e <predict on minus strand only>

  -n <name for the gene [default snap]>

  -i <FASTA file of proteins>
  
  -t <FASTA file of transcripts>

  -l <treat lowercase as N>

  -h Show this usage information

EOF
    exit 0
}

plusmode=0
minusmode=0

while getopts ":m:f:g:ha:pen:i:t:" opt; do
  case $opt in
    m)
     hmmfile=$OPTARG
      ;;
    f)
     fastafile=$OPTARG
      ;;
    h)
     usage
     exit 1
      ;;    
    g)
     gfffile=$OPTARG
      ;;
    a)
     acefile=$OPTARG
      ;;  
    p)
     plusmode=$OPTARG
      ;;
    e)
     minusmode=$OPTARG
      ;;
    n)
     name=$OPTARG
      ;;  
    i)
     aminofile=$OPTARG
      ;;
    t)
     transcripts=$OPTARG
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


if [ ! -z "$gfffile" ] && [ "$plusmode" == 0 ] && [ "$minusmode" == 0 ]; then
     snap -quiet -lcmask -name "$name" "$hmmfile" "$fastafile" -gff -aa "$aminofile" -tx "$transcripts" > "$gfffile"

elif [ ! -z "$gfffile" ] && [ "$plusmode" != 0 ] && [ "$minusmode" == 0 ]; then
     snap -quiet -lcmask -plus -name "$name" "$hmmfile" "$fastafile" -gff -aa "$aminofile" -tx "$transcripts" > "$gfffile"

elif [ ! -z "$gfffile" ] && [ "$plusmode" == 0 ] && [ "$minusmode" != 0 ]; then
     snap -quiet -lcmask -minus -name "$name" "$hmmfile" "$fastafile" -gff -aa "$aminofile" -tx "$transcripts" > "$gfffile"

elif [ ! -z $acefile ]; then
     snap -quiet -lcmask -name "$name" "$hmmfile" "$fastafile" -ace -aa "$aminofile" -tx "$transcripts" > "$acefile"

fi
