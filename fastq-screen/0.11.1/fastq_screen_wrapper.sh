#!/bin/bash

set -e
set -x

usage() {
      echo ""
      echo "Usage : sh $0 -g <reference_genome> [-I <Index_folder>] -U <reads> -A <aligner> -t <number> -s <subset> -o <output>"
      echo "Example: sh fastq_screen_wrapper.sh -I new -I new2 -U fqs_test_dataset.fastq.gz -A bowtie2 -t 6 -s 100000 -o final"
      echo ""

cat <<'EOF'
  
  ###### Command line options ##########

  -g                    Specify the reference genome fasta file. If there are more than one specify 
                        them like -g genome1 -g genome2...
  -I                    Specify the reference genomeindex folder if you already have a reference 
                        genome indexed with bowtie, bowtie2 and bwa. If there are more than one specify 
                        them like -I index1 -I index2...
  -U                    Specify ead files. Both compressed and uncompressed works and more than one 
                        read files can be specified
  -A                    Specify the aligner to use for the mapping. Valid arguments are 'bowtie', bowtie2' 
                        or 'bwa'. Only one specifier can be specified at a time
  -t                    Number of threads (Default is 1)                 
  -i                    Assume that the quality values are in encoded in Illumina v1.3 format. Defaults 
                        to Sanger format if this flag is not specified.
  -s                    By Default FastQ Screen runs with this parameter set to 100000. To process
                        an entire dataset however, adjust --subset to 0
  -o                    Specify a directory in which to save output files.If no directory is specified then 
                        output files are saved into the same directory as the input file.
EOF
    exit 0
}

illumina1_3=0

while getopts ":hg:U:A:t:io:I:s:" opt; do
  case $opt in
    g)
    referencegenome+=("$OPTARG")
     ;;
    U)
    reads+=("$OPTARG")
     ;;
    A)
    aligner=$OPTARG 
     ;;
    t)
    threads=$OPTARG
     ;;
    o)
    output=$OPTARG
     ;;
    i)
    illumina1_3=$OPTARG
     ;;
    I)
    index+=("$OPTARG")
     ;;
    s)
     subset=$OPTARG
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

if [ ! -z "$referencegenome" ] && [ ! -z "$index" ]; then # Reference Genomes and Indexes both

    for f in "${referencegenome[@]}"; do
        filename=$(basename "$f")
        extension="${filename##*.}"
        filename="${filename%.*}"
        if [[ "$extension" =~ "gz" ]]; then
          gunzip -c $f > $filename
          extension="${filename##*.}"
          prefix="${filename%.*}"
          mkdir $prefix
          if [[ "$aligner" == "bowtie" ]]; then
            bowtie-build $filename $prefix
            mv $prefix*ebwt $prefix
            echo "DATABASE        $prefix $prefix/$prefix BOWTIE" >> fastq_screen.conf
          elif [[ "$aligner" == "bowtie2" ]]; then
            bowtie2-build $filename $prefix
            mv $prefix*bt2 $prefix
            echo "DATABASE        $prefix $prefix/$prefix BOWTIE2" >> fastq_screen.conf
          elif [[ "$aligner" == "bwa" ]]; then
            bwa index $filename -p $prefix
            mv $prefix*bwt $prefix*pac $prefix*ann $prefix*amb $prefix*sa $prefix
            echo "DATABASE        $prefix $prefix/$prefix BWA" >> fastq_screen.conf
          fi
        elif [[ "$extension" != "gz" ]]; then
          mkdir $filename
          if [[ "$aligner" == "bowtie" ]]; then
            bowtie-build $f $filename
            mv $filename*ebwt $filename
            echo "DATABASE        $filename $filename/$filename BOWTIE" >> fastq_screen.conf
          elif [[ "$aligner" == "bowtie2" ]]; then
            bowtie2-build $f $filename
            mv $filename*bt2 $filename
            echo "DATABASE        $filename $filename/$filename BOWTIE2" >> fastq_screen.conf
          elif [[ "$aligner" == "bwa" ]]; then
            bwa index $f -p $filename
            mv $filename*bwt $filename*pac $filename*ann $filename*amb $filename*sa $filename
            echo "DATABASE        $filename $filename/$filename BWA" >> fastq_screen.conf
          fi
        fi
    done

      if [[ "$aligner" == "bowtie" ]]; then
        for x in "${index[@]}"; do
          fil=$(ls $x/*ebwt | head -n 1)
          filename="${fil%.*.*}"
          ind=$(echo $filename | cut -d "/" -f 1 )
          px=$(echo $filename | cut -d "/" -f 2 )
          echo "DATABASE        $px $ind/$px BOWTIE" >> fastq_screen.conf
        done
      elif [[ "$aligner" == "bowtie2" ]]; then
        for x in "${index[@]}"; do
          fil=$(ls $x/*bt2 | head -n 1)
          filename="${fil%.*.*}"
          ind=$(echo $filename | cut -d "/" -f 1 )
          px=$(echo $filename | cut -d "/" -f 2 )
          echo "DATABASE        $px $ind/$px BOWTIE2" >> fastq_screen.conf
        done
      elif [[ "$aligner" == "bwa" ]]; then
        for x in "${index[@]}"; do
          fil=$(ls $x/*bwt | head -n 1)
          filename="${fil%.*.*}"
          ind=$(echo $filename | cut -d "/" -f 1 )
          px=$(echo $filename | cut -d "/" -f 2 )
          echo "DATABASE        $px $ind/$px BWA" >> fastq_screen.conf
        done
      fi

elif [ ! -z "$referencegenome" ] && [ -z "$index" ]; then # Reference Genomes only


    for f in "${referencegenome[@]}"; do
        filename=$(basename "$f")
        extension="${filename##*.}"
        filename="${filename%.*}"
        if [[ "$extension" =~ "gz" ]]; then
          gunzip -c $f > $filename
          extension="${filename##*.}"
          prefix="${filename%.*}"
          mkdir $prefix
          if [[ "$aligner" == "bowtie" ]]; then
            bowtie-build $filename $prefix
            mv $prefix*ebwt $prefix
            echo "DATABASE        $prefix $prefix/$prefix BOWTIE" >> fastq_screen.conf
          elif [[ "$aligner" == "bowtie2" ]]; then
            bowtie2-build $filename $prefix
            mv $prefix*bt2 $prefix
            echo "DATABASE        $prefix $prefix/$prefix BOWTIE2" >> fastq_screen.conf
          elif [[ "$aligner" == "bwa" ]]; then
            bwa index $filename -p $prefix
            mv $prefix*bwt $prefix*pac $prefix*ann $prefix*amb $prefix*sa $prefix
            echo "DATABASE        $prefix $prefix/$prefix BWA" >> fastq_screen.conf
          fi
        elif [[ "$extension" != "gz" ]]; then
          mkdir $filename
          if [[ "$aligner" == "bowtie" ]]; then
            bowtie-build $f $filename
            mv $filename*ebwt $filename
            echo "DATABASE        $filename $filename/$filename BOWTIE" >> fastq_screen.conf
          elif [[ "$aligner" == "bowtie2" ]]; then
            bowtie2-build $f $filename
            mv $filename*bt2 $filename
            echo "DATABASE        $filename $filename/$filename BOWTIE2" >> fastq_screen.conf
          elif [[ "$aligner" == "bwa" ]]; then
            bwa index $f -p $filename
            mv $filename*bwt $filename*pac $filename*ann $filename*amb $filename*sa $filename
            echo "DATABASE        $filename $filename/$filename BWA" >> fastq_screen.conf
          fi
        fi
    done

elif [  -z "$referencegenome" ] && [ ! -z "$index" ]; then # Indexes only

      if [[ "$aligner" == "bowtie" ]]; then
        for x in "${index[@]}"; do
          fil=$(ls $x/*ebwt | head -n 1)
          filename="${fil%.*.*}"
          ind=$(echo $filename | cut -d "/" -f 1 )
          px=$(echo $filename | cut -d "/" -f 2 )
          echo "DATABASE        $px $ind/$px BOWTIE" >> fastq_screen.conf
        done
      elif [[ "$aligner" == "bowtie2" ]]; then
        for x in "${index[@]}"; do
          fil=$(ls $x/*bt2 | head -n 1)
          filename="${fil%.*.*}"
          ind=$(echo $filename | cut -d "/" -f 1 )
          px=$(echo $filename | cut -d "/" -f 2 )
          echo "DATABASE        $px $ind/$px BOWTIE2" >> fastq_screen.conf
        done
      elif [[ "$aligner" == "bwa" ]]; then
        for x in "${index[@]}"; do
          fil=$(ls $x/*bwt | head -n 1)
          filename="${fil%.*.*}"
          ind=$(echo $filename | cut -d "/" -f 1 )
          px=$(echo $filename | cut -d "/" -f 2 )
          echo "DATABASE        $px $ind/$px BWA" >> fastq_screen.conf
        done
      fi
fi

# Fastq_screen run

if [ "$illumina1_3" == 0 ]; then
  fastq_screen --conf fastq_screen.conf $reads --aligner $aligner --threads $threads --outdir $output --subset $subset
elif [ "$illumina1_3" != 0 ]; then
  fastq_screen --conf fastq_screen.conf $reads --aligner $aligner --threads $threads --outdir $output --illumina1_3 --subset $subset
fi