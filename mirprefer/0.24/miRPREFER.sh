#!/bin/bash

usage() {
       echo ""
      echo "Usage : sh $0 -f FASTA_FILE -a ALIGNMENT_FILE -G GFF_FILE_EXCLUDE -I GFF_FILE_INCLUDE -p PRECURSOR_LEN -r READS_DEPTH_CUTOFF -c NUM_OF_CORE -o OUTFOLDER -i NAME_PREFIX -g MAX_GAP -m MIN_MATURE_LEN -x MAX_MATURE_LEN -s ALLOW_NO_STAR_EXPRESSION -t ALLOW_3NT_OVERHANG"
      echo ""

cat <<'EOF'

  -f </path/to/fasta file>
  -a </path/to/alignment file>
  -G GFF_FILE_EXCLUDE # This option is mutual exclusive with 'GFF_FILE_INCLUDE', Thus only one of them can be use.
  -I GFF_FILE_INCLUDE # This option is mutual exclusive with 'GFF_FILE_EXCLUDE'. Thus, only one of them can be used.
  -p PRECURSOR_LEN
  -r READS_DEPTH_CUTOFF
  -c NUM_OF_CORE
  -o OUTFOLDER
  -i NAME_PREFIX
  -g MAX_GAP
  -m MIN_MATURE_LEN
  -x MAX_MATURE_LEN
  -s ALLOW_NO_STAR_EXPRESSION
  -t ALLOW_3NT_OVERHANG

EOF
    exit 0
}

while getopts ":f:a:G:I:p:hr:c:o:i:g:m:x:s:t:" opt; do
  case $opt in
    f)
     fasta=$OPTARG
      ;;
    a)
     alignmentfile+=("$OPTARG");;
    h)
     usage
     exit 1
      ;;    
    G)
     gfffileex=$OPTARG
      ;;
    I)
     gfffilein=$OPTARG
      ;;
    p)
     plength=$OPTARG # This option is specific to CyVerse Discovery Environment
      ;;
    r)
     readdpthcut=$OPTARG
      ;;
    c)
     cores=$OPTARG
      ;;  
    o)
     output=$OPTARG
      ;;
    i)
     prefix=$OPTARG
      ;;
    g)
     maxgap=$OPTARG
      ;;
    m)
     minmatl=$OPTARG
      ;;
    x)
     maxmatl=$OPTARG
      ;;
    s)
     nostar=$OPTARG
      ;;
    t)
     overhang=$OPTARG
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

new=$(printf '%s,' "${alignmentfile[@]}")
var1=$(sed 's/,$//' <<< "$new")

echo "PIPELINE_PATH = /miR-PREFeR-0.24" >> test.config
echo "FASTA_FILE =  $fasta" >> test.config
echo "ALIGNMENT_FILE = $var1" >> test.config
if [[ ! -z $gfffileex ]]; then
	echo "GFF_FILE_EXCLUDE = $gfffileex" >> test.config
elif [[ ! -z $gfffilein ]]; then
	echo "GFF_FILE_INCLUDE = $gfffilein" >> test.config
fi
echo "PRECURSOR_LEN = $plength" >> test.config
echo "READS_DEPTH_CUTOFF = $readdpthcut" >> test.config
echo "NUM_OF_CORE = $cores" >> test.config
echo "OUTFOLDER = $output" >> test.config
echo "NAME_PREFIX = $prefix" >> test.config
echo "MAX_GAP = $maxgap" >> test.config
echo "MIN_MATURE_LEN = $minmatl" >> test.config
echo "MAX_MATURE_LEN = $maxmatl" >> test.config
echo "ALLOW_NO_STAR_EXPRESSION = $nostar" >> test.config
echo "ALLOW_3NT_OVERHANG = $overhang" >> test.config
echo "CHECKPOINT_SIZE = 3000" >> test.config

# miRPREFER python script
python /miR-PREFeR-0.24/miR_PREFeR.py -L -k pipeline test.config
