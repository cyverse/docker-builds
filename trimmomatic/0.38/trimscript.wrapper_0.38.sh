#!/bin/bash

set -x
set -e

quality_64=0

while getopts ":h:1:2:U:p:r:a:QO:s:c:t:" opt; do
  case $opt in
    1)
    left_reads+=("$OPTARG") # Left reads
     ;;
    2)
    right_reads=("$OPTARG") # Right reads
     ;;
    U)
    single_reads+=("$OPTARG") # Single end reads;;
     ;;
    O)
    outdir=$OPTARG # Output directory
     ;;
    p)
    num_threads=$OPTARG # Number of threads
     ;;
    r)
    trimfile=$OPTARG # Trimmer setting file
     ;;
    a)
    adapter_fle=$OPTARG # Adpater file
     ;;
    Q)
    quality_64=$OPTARG # Phread 64
     ;;
    s)
    seed=$OPTARG # Seed mismatches. Recommended value is 2 (Single end)
     ;;
    c)
    clip=$OPTARG # palindromeClipThreshold. Values around 30 or more are recommended 
     ;;
    t)
    simple=$OPTARG # simpleClipThreshold. Values between 7 and 15 are recommended but this depends on the length of the reads (Single end). 10 may be
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

cat $trimfile | sed 's/ILLUMINACLIP:/ILLUMINACLIP:\/staging\//' > trimfile.txt

mkdir $outdir

cuts=`cat trimfile.txt | grep -v \#`

rm trimfile.txt

if [ "$quality_64" == 0 ]; then

  if [ ! -z "$left_reads" ] && [ ! -z "$right_reads" ]; then
      for f in "${left_reads[@]}"; do
        extension=$(echo "$f" | sed -r 's/.*(fq|fq.gz|fastq|fastq.gz)$/\1/')
        if [[ "$extension" =~ "fq.gz" ]]; then
          filename=$(basename "$f" ".fq.gz")
      	  filename2=${filename/_R1/_R2}
          trimmomatic PE -threads $num_threads ${filename}.fq.gz ${filename2}.fq.gz 'trmPr_'${filename}.fq.gz 'trmS_'${filename}.fq.gz 'trmPr_'${filename2}.fq.gz 'trmS_'${filename2}.fq.gz -phred64 ILLUMINACLIP:"$adapter_fle":"$seed":"$clip":"$simple" $cuts
          mv 'trmPr_'${filename}.fq.gz 'trmS_'${filename}.fq.gz 'trmPr_'${filename2}.fq.gz 'trmS_'${filename2}.fq.gz $outdir
        elif [[ "$extension" =~ "fastq.gz" ]]; then
          filename=$(basename "$f" ".fastq.gz")
          filename2=${filename/_R1/_R2}
          trimmomatic PE -threads $num_threads ${filename}.fastq.gz ${filename2}.fastq.gz 'trmPr_'${filename}.fastq.gz 'trmS_'${filename}.fastq.gz 'trmPr_'${filename2}.fastq.gz 'trmS_'${filename2}.fastq.gz -phred64 ILLUMINACLIP:"$adapter_fle":"$seed":"$clip":"$simple" $cuts
          mv 'trmPr_'${filename}.fastq.gz 'trmS_'${filename}.fastq.gz 'trmPr_'${filename2}.fastq.gz 'trmS_'${filename2}.fastq.gz $outdir
        elif [[ "$extension" =~ "fq" ]]; then
          filename=$(basename "$f" ".fq")
          filename2=${filename/_R1/_R2}
          trimmomatic PE -threads $num_threads ${filename}.fq ${filename2}.fq 'trmPr_'${filename}.fq 'trmS_'${filename}.fq 'trmPr_'${filename2}.fq 'trmS_'${filename2}.fq -phred64 ILLUMINACLIP:"$adapter_fle":"$seed":"$clip":"$simple" $cuts
          mv 'trmPr_'${filename}.fq 'trmS_'${filename}.fq 'trmPr_'${filename2}.fq 'trmS_'${filename2}.fq $outdir
        elif [[ "$extension" =~ "fastq" ]]; then
          filename=$(basename "$f" ".fastq")
          filename2=${filename/_R1/_R2}
          trimmomatic PE -threads $num_threads ${filename}.fastq ${filename2}.fastq 'trmPr_'${filename}.fastq 'trmS_'${filename}.fastq 'trmPr_'${filename2}.fastq 'trmS_'${filename2}.fastq -phred64 ILLUMINACLIP:"$adapter_fle":"$seed":"$clip":"$simple" $cuts
          mv 'trmPr_'${filename}.fastq 'trmS_'${filename}.fastq 'trmPr_'${filename2}.fastq 'trmS_'${filename2}.fastq $outdir
        elif [ "$extension" != "fastq" ] || [ "$extension" != "fq" ] || [ "$extension" != "fastq.gz" ] || [ "$extension" != "fq.gz" ]; then
          echo "The extension" "$extension" "is not supported. Only .fq, .fq.gz, .fastq, .fastq.gz are only supported" 1>&2        
          exit 64
        fi 
      done

  # Single end reads

  elif [ ! -z "$single_reads" ]; then
      for f in "${single_reads[@]}"; do
        trimmomatic SE -threads $num_threads $f 'trimS_'$f -phred64 ILLUMINACLIP:"$adapter_fle":"$seed":"$clip":"$simple" $cuts
        mv trimS_* $outdir
      done
  fi

elif [ "$quality_64" != 0 ]; then

    if [ ! -z "$left_reads" ] && [ ! -z "$right_reads" ]; then
      for f in "${left_reads[@]}"; do
        extension=$(echo "$f" | sed -r 's/.*(fq|fq.gz|fastq|fastq.gz)$/\1/')
        if [[ "$extension" =~ "fq.gz" ]]; then
          filename=$(basename "$f" ".fq.gz")
          filename2=${filename/_R1/_R2}
          trimmomatic PE -threads $num_threads ${filename}.fq.gz ${filename2}.fq.gz 'trmPr_'${filename}.fq.gz 'trmS_'${filename}.fq.gz 'trmPr_'${filename2}.fq.gz 'trmS_'${filename2}.fq.gz -phred33 ILLUMINACLIP:"$adapter_fle":"$seed":"$clip":"$simple" $cuts
          mv 'trmPr_'${filename}.fq.gz 'trmS_'${filename}.fq.gz 'trmPr_'${filename2}.fq.gz 'trmS_'${filename2}.fq.gz $outdir
        elif [[ "$extension" =~ "fastq.gz" ]]; then
          filename=$(basename "$f" ".fastq.gz")
          filename2=${filename/_R1/_R2}
          trimmomatic PE -threads $num_threads ${filename}.fastq.gz ${filename2}.fastq.gz 'trmPr_'${filename}.fastq.gz 'trmS_'${filename}.fastq.gz 'trmPr_'${filename2}.fastq.gz 'trmS_'${filename2}.fastq.gz -phred33 ILLUMINACLIP:"$adapter_fle":"$seed":"$clip":"$simple" $cuts
          mv 'trmPr_'${filename}.fastq.gz 'trmS_'${filename}.fastq.gz 'trmPr_'${filename2}.fastq.gz 'trmS_'${filename2}.fastq.gz $outdir
        elif [[ "$extension" =~ "fq" ]]; then
          filename=$(basename "$f" ".fq")
          filename2=${filename/_R1/_R2}
          trimmomatic PE -threads $num_threads ${filename}.fq ${filename2}.fq 'trmPr_'${filename}.fq 'trmS_'${filename}.fq 'trmPr_'${filename2}.fq 'trmS_'${filename2}.fq -phred33 ILLUMINACLIP:"$adapter_fle":"$seed":"$clip":"$simple" $cuts
          mv 'trmPr_'${filename}.fq 'trmS_'${filename}.fq 'trmPr_'${filename2}.fq 'trmS_'${filename2}.fq $outdir
        elif [[ "$extension" =~ "fastq" ]]; then
          filename=$(basename "$f" ".fastq")
          filename2=${filename/_R1/_R2}
          trimmomatic PE -threads $num_threads ${filename}.fastq ${filename2}.fastq 'trmPr_'${filename}.fastq 'trmS_'${filename}.fastq 'trmPr_'${filename2}.fastq 'trmS_'${filename2}.fastq -phred33 ILLUMINACLIP:"$adapter_fle":"$seed":"$clip":"$simple" $cuts
          mv 'trmPr_'${filename}.fastq 'trmS_'${filename}.fastq 'trmPr_'${filename2}.fastq 'trmS_'${filename2}.fastq $outdir
        elif [ "$extension" != "fastq" ] || [ "$extension" != "fq" ] || [ "$extension" != "fastq.gz" ] || [ "$extension" != "fq.gz" ]; then
          echo "The extension" "$extension" "is not supported. Only .fq, .fq.gz, .fastq, .fastq.gz are only supported" 1>&2        
          exit 64
        fi 
      done

  # Single end reads

  elif [ ! -z "$single_reads" ]; then
      for f in "${single_reads[@]}"; do
        trimmomatic SE -threads $num_threads $f 'trimS_'$f -phred33 ILLUMINACLIP:"$adapter_fle":"$seed":"$clip":"$simple" $cuts
        mv trimS_* $outdir
      done
  fi

fi
