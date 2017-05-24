#!/bin/bash

usage() {
      echo ""
      echo "Usage : sh $0 -f SAMTYPE -r ORDER -s STRANDED -a MINAQUAL -t FEATURETYPE -i IDATTR -m MODE -o <sam_bam out> -n <sam/bam files> -g <annotation file gff|gtf> -c <counts folder>" 
      echo ""

cat <<'EOF'
  -f type of <alignment_file> data, either 'sam' or 'bam'
                        (default: sam)  

  -r 'pos' or 'name'. Sorting order of <alignment_file>
                        (default: name). Paired-end sequencing data must be
                        sorted either by position or by read name, and the
                        sorting order must be specified. Ignored for single-
                        end data.

  -s whether the data is from a strand-specific assay.
                        Specify 'yes', 'no', or 'reverse' (default: yes).
                        'reverse' means 'yes' with reversed strand
                        interpretation

  -a skip all reads with alignment quality lower than the
                        given minimum value (default: 10)


  -t feature type (3rd column in GFF file) to be used, all
                        features of other type are ignored (default, suitable
                        for Ensembl GTF files: exon)

  -i GFF attribute to be used as feature ID (default,
                        suitable for Ensembl GTF files: gene_id)

  -m mode to handle reads overlapping more than one feature
                        (choices: union, intersection-strict, intersection-
                        nonempty; default: union)

  -n folder containing sam or bam input files

  -g gff or gtf file

  -c output folder containing the final counts files
                        
EOF
    exit 0
}

while getopts ":hf:r:a:s:t:i:m:n:g:o:c:" opt; do
  case $opt in
    f)
    types=$OPTARG 
     ;;
    r)
    pos=$OPTARG 
     ;;
    a)
    qual=$OPTARG
     ;;
    s)
    strand=$OPTARG 
     ;;
    t)
    feature=$OPTARG 
     ;;
    i)
    attribute=$OPTARG 
     ;;
    m)
    mode=$OPTARG 
     ;;
    o)
    sam_out=$OPTARG
     ;; 
    n)
     multi+=("$OPTARG")
     ;;
    g)
    input_gff=$OPTARG
     ;;
    c)
    counts_folder=$OPTARG
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

shift $((OPTIND -1))

mkdir "$counts_folder"

function no_out_sam()
{
    filename=$(basename "$file")
    extension="${file##*.}"
              
    if [ $extension = "sam" ];
      then
        filename="${filename%.*}"
        feat_out=$filename.out
        htseq-count -f $types -r "$pos" -m "$mode" -i "$attribute" -s "$strand" -t "$feature" -a "$qual" "$file" "$input_gff" > "$feat_out"
        mv "$feat_out" "$counts_folder"
    elif [ $extension = "bam" ];
      then
        filename="${filename%.*}"
        feat_out=$filename.out
        htseq-count -f $types -r "$pos" -m "$mode" -i "$attribute" -s "$strand" -t "$feature" -a "$qual" "$file" "$input_gff" > "$feat_out"
        mv "$feat_out" "$counts_folder"
    fi    
}


for file in "${multi[@]}"; do 
    no_out_sam &
done
wait

