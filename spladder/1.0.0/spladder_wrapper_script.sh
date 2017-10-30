#!/bin/bash

set -e
set -x

while getopts ":hp:b:o:a:t:c:P:M:" opt; do
  case $opt in
        p)
        threads=$OPTARG
        ;;
        b)
        bams=$OPTARG
        ;;
        o)
        outdir=$OPTARG
        ;;
        a)
        anno=$OPTARG
        ;;
        t)
        event+=("$OPTARG")
        ;;
        c)
        confidence=$OPTARG
        ;;
        P)
        primary=$OPTARG
        ;;
        M)
        merge=$OPTARG
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

### set dirs and variables
bams=$(ls $bams/*.bam | tr "\n" "," | sed 's/,$//')

if [[ "$event" =~ "all" ]]; then
    python /spladder-1.0.0/python/spladder.py -b $bams -o $outdir -a $anno -v y -M $merge -t exon_skip,intron_retention,alt_3prime,alt_5prime,mutex_exons,mult_exon_skip -c $confidence -p $threads -P $primary
else
    python /spladder-1.0.0/python/spladder.py -b $bams -o $outdir -a $anno -v y -M $merge -t $event -c $confidence -p $threads -P $primary
fi


