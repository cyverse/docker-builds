#!/bin/bash

set -x
set -e

while getopts ":hg:1:2:r:n:p:t:a:b:c:d:e:f:i:j:k:l:s:m:" opt; do
  case $opt in
    g)
    referencegenome=$OPTARG # Reference genome file
     ;;
    r)
    referenceannotation=$OPTARG # Reference annotation
     ;;
    1)
    left_reads+=("$OPTARG") # Left reads
     ;;
    2)
    right_reads=("$OPTARG") # Right reads
     ;;
    n)
    threads=$OPTARG
     ;;
    p)
    repeat=$OPTARG
     ;;
    t)
    strand=$OPTARG
     ;;
    a)
    bam+=("$OPTARG")
     ;;
    c)
    first=$OPTARG
     ;;
    d)
    cluster_distance=$OPTARG
     ;;
    s)
    cluster_size_alu_ad1=$OPTARG
     ;;
    b)
    cluster_size_alu_ad2=$OPTARG
     ;;
    e)
    cluster_size_nalurp=$OPTARG
     ;;
    f)
    cluster_size_nrp=$OPTARG
     ;;
    i)
    cluster_size_rg=$OPTARG
     ;;
    j)
    cluster_size_hp=$OPTARG
     ;;
    k)
    cluster_size_alu_hp=$OPTARG
     ;;
    l)
    cluster_size_nalurp_hp=$OPTARG
     ;;
    m)
    cluster_size_nrp_hp=$OPTARG
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

if [ ! -z "$referenceannotation" ]; then
    sprint prepare -t $referenceannotation $referencegenome bwa
else
    sprint prepare $referencegenome bwa
fi

if [ ! -z "$left_reads" ] && [ ! -z "$right_reads" ]; then
    
    for f in "${left_reads[@]}"; do
        extension=$(echo "$f" | sed -r 's/.*(fq|fq.gz|fastq|fastq.gz)$/\1/')
        if [[ "$extension" =~ "fastq.gz" ]]; then
            filename=$(basename "$f" ".fastq.gz")
            filename2=${filename/R1/R2}
            filename3=$(echo $filename | sed 's/R1_//')
            bwa mem -t $threads $referencegenome ${filename}.fastq.gz ${filename2}.fastq.gz > $filename3.sam
            samtools view -bS $filename3.sam > $filename3.bam
            samtools sort $filename3.bam $filename3.sorted
            if [ ! -z "$repeat" ]; then
                sprint main -c $first -cd $cluster_distance -csad1 $cluster_size_alu_ad1 -csad2 $cluster_size_alu_ad2 -csnar $cluster_size_nalurp -csnr $cluster_size_nrp -csahp $cluster_size_hp -csnarhp $cluster_size_alu_hp -csnrhp $cluster_size_nalurp_hp -rp $repeat -ss $strand -p $threads -1 ${filename}.fastq.gz -2 ${filename2}.fastq.gz $referencegenome $filename3"_SPRINT_output" bwa samtools 
            else
                sprint main -c $first -cd $cluster_distance -csrg $cluster_size_rg -cshp $cluster_size_nrp_hp -p $threads -ss $strand -1 ${filename}.fastq.gz -2 ${filename2}.fastq.gz $referencegenome $filename3"_SPRINT_output" bwa samtools
            fi
            #rm $filename3.sam
        elif [[ "$extension" =~ "fq.gz" ]]; then
            filename=$(basename "$f" ".fq.gz")
            filename2=${filename/R1/R2}
            filename3=$(echo $filename | sed 's/R1_//')
            bwa mem $referencegenome ${filename}.fq.gz ${filename2}.fq.gz > $filename3.sam
            samtools view -bS $filename3.sam > $filename3.bam
            samtools sort $filename3.bam $filename3.sorted
            if [ ! -z "$repeat" ]; then
                sprint main -c $first -cd $cluster_distance -csad1 $cluster_size_alu_ad1 -csad2 $cluster_size_alu_ad2 -csnar $cluster_size_nalurp -csnr $cluster_size_nrp -csahp $cluster_size_hp -csnarhp $cluster_size_alu_hp -csnrhp $cluster_size_nalurp_hp -rp $repeat -ss $strand -p $threads -1 ${filename}.fq.gz -2 ${filename2}.fq.gz $referencegenome $filename3"_SPRINT_output" bwa samtools
            else
                sprint main -c $first -cd $cluster_distance -csrg $cluster_size_rg -cshp $cluster_size_nrp_hp -p $threads -ss $strand -1 ${filename}.fq.gz -2 ${filename2}.fq.gz $referencegenome $filename3"_SPRINT_output" bwa samtools
            fi
            rm $filename3.sam
        elif [[ "$extension" =~ "fastq" ]]; then
            filename=$(basename "$f" ".fastq")
            filename2=${filename/R1/R2}
            filename3=$(echo $filename | sed 's/R1_//')
            bwa mem $referencegenome ${filename}.fastq ${filename2}.fastq > $filename3.sam
            samtools view -bS $filename3.sam > $filename3.bam
            samtools sort $filename3.bam $filename3.sorted
            if [ ! -z "$repeat" ]; then
                sprint main -c $first -cd $cluster_distance -csad1 $cluster_size_alu_ad1 -csad2 $cluster_size_alu_ad2 -csnar $cluster_size_nalurp -csnr $cluster_size_nrp -csahp $cluster_size_hp -csnarhp $cluster_size_alu_hp -csnrhp $cluster_size_nalurp_hp -rp $repeat -ss $strand -p $threads -1 ${filename}.fastq -2 ${filename2}.fastq $referencegenome $filename3"_SPRINT_output" bwa samtools
            else
                sprint main -c $first -cd $cluster_distance -csrg $cluster_size_rg -cshp $cluster_size_nrp_hp -p $threads -ss $strand -1 ${filename}.fastq -2 ${filename2}.fastq $referencegenome $filename3"_SPRINT_output" bwa samtools
            fi
            rm $filename3.sam
        elif [[ "$extension" =~ "fq" ]]; then
            filename=$(basename "$f" ".fq")
            filename2=${filename/R1/R2}
            filename3=$(echo $filename | sed 's/R1_//')
            bwa mem $referencegenome ${filename}.fq ${filename2}.fq > $filename3.sam
            samtools view -bS $filename3.sam > $filename3.bam
            samtools sort $filename3.bam $filename3.sorted
            if [ ! -z "$repeat" ]; then
                sprint main -c $first -cd $cluster_distance -csad1 $cluster_size_alu_ad1 -csad2 $cluster_size_alu_ad2 -csnar $cluster_size_nalurp -csnr $cluster_size_nrp -csahp $cluster_size_hp -csnarhp $cluster_size_alu_hp -csnrhp $cluster_size_nalurp_hp -rp $repeat -ss $strand -p $threads -1 ${filename}.fq -2 ${filename2}.fq $referencegenome $filename3"_SPRINT_output" bwa samtools
            else
                sprint main -c $first -cd $cluster_distance -csrg $cluster_size_rg -cshp $cluster_size_nrp_hp -p $threads -ss $strand -1 ${filename}.fq -2 ${filename2}.fq $referencegenome $filename3"_SPRINT_output" bwa samtools
            fi
            rm $filename3.sam
        fi
    done

elif [ ! -z "$left_reads" ] && [ -z "$right_reads" ]; then

    for f in "${left_reads[@]}"; do
        extension=$(echo "$f" | sed -r 's/.*(fq|fq.gz|fastq|fastq.gz)$/\1/')
        if [[ "$extension" =~ "fastq.gz" ]]; then
            filename=$(basename "$f" ".fastq.gz")
            filename2=$(echo $filename | sed 's/R1_//')
            bwa mem -t $threads $referencegenome ${filename}.fastq.gz > $filename2.sam
            samtools view -bS $filename2.sam > $filename2.bam
            samtools sort $filename2.bam $filename2.sorted
            if [ ! -z "$repeat" ]; then
                sprint main -c $first -cd $cluster_distance -csad1 $cluster_size_alu_ad1 -csad2 $cluster_size_alu_ad2 -csnar $cluster_size_nalurp -csnr $cluster_size_nrp -csahp $cluster_size_hp -csnarhp $cluster_size_alu_hp -csnrhp $cluster_size_nalurp_hp -rp $repeat -ss $strand -p $threads -1 ${filename}.fastq.gz $referencegenome $filename2"_SPRINT_output" bwa samtools
            else
                sprint main -c $first -cd $cluster_distance -csrg $cluster_size_rg -cshp $cluster_size_nrp_hp -p $threads -ss $strand -1 ${filename}.fastq.gz $referencegenome $filename2"_SPRINT_output" bwa samtools
            fi
            rm $filename2.sam
        elif [[ "$extension" =~ "fq.gz" ]]; then
            filename=$(basename "$f" ".fq.gz")
            filename2=$(echo $filename | sed 's/R1_//')
            bwa mem $referencegenome ${filename}.fq.gz > $filename2.sam
            samtools view -bS $filename2.sam > $filename2.bam
            samtools sort $filename2.bam $filename2.sorted
            if [ ! -z "$repeat" ]; then
                sprint main -c $first -cd $cluster_distance -csad1 $cluster_size_alu_ad1 -csad2 $cluster_size_alu_ad2 -csnar $cluster_size_nalurp -csnr $cluster_size_nrp -csahp $cluster_size_hp -csnarhp $cluster_size_alu_hp -csnrhp $cluster_size_nalurp_hp -rp $repeat -ss $strand -p $threads -1 ${filename}.fq.gz $referencegenome $filename2"_SPRINT_output" bwa samtools
            else
                sprint main -c $first -cd $cluster_distance -csrg $cluster_size_rg -cshp $cluster_size_nrp_hp -p $threads -ss $strand -1 ${filename}.fq.gz $referencegenome $filename2"_SPRINT_output" bwa samtools
            fi
            rm $filename2.sam
        elif [[ "$extension" =~ "fastq" ]]; then
            filename=$(basename "$f" ".fastq")
            filename2=$(echo $filename | sed 's/R1_//')
            bwa mem $referencegenome ${filename}.fastq > $filename2.sam
            samtools view -bS $filename2.sam > $filename2.bam
            samtools sort $filename2.bam $filename2.sorted
            if [ ! -z "$repeat" ]; then
                sprint main -c $first -cd $cluster_distance -csad1 $cluster_size_alu_ad1 -csad2 $cluster_size_alu_ad2 -csnar $cluster_size_nalurp -csnr $cluster_size_nrp -csahp $cluster_size_hp -csnarhp $cluster_size_alu_hp -csnrhp $cluster_size_nalurp_hp -rp $repeat -ss $strand -p $threads -1 ${filename}.fastq $referencegenome $filename2"_SPRINT_output" bwa samtools
            else
                sprint main -c $first -cd $cluster_distance -csrg $cluster_size_rg -cshp $cluster_size_nrp_hp -p $threads -ss $strand -1 ${filename}.fastq $referencegenome $filename2"_SPRINT_output" bwa samtools
            fi
            rm $filename2.sam
        elif [[ "$extension" =~ "fq" ]]; then
            filename=$(basename "$f" ".fq")
            filename2=$(echo $filename | sed 's/R1_//')
            bwa mem $referencegenome ${filename}.fq > $filename2.sam
            samtools view -bS $filename2.sam > $filename2.bam
            samtools sort $filename2.bam $filename2.sorted
            if [ ! -z "$repeat" ]; then
                sprint main -c $first -cd $cluster_distance -csad1 $cluster_size_alu_ad1 -csad2 $cluster_size_alu_ad2 -csnar $cluster_size_nalurp -csnr $cluster_size_nrp -csahp $cluster_size_hp -csnarhp $cluster_size_alu_hp -csnrhp $cluster_size_nalurp_hp -rp $repeat -ss $strand -p $threads -1 ${filename}.fq $referencegenome $filename2"_SPRINT_output" bwa samtools
            else
                sprint main -c $first -cd $cluster_distance -csrg $cluster_size_rg -cshp $cluster_size_nrp_hp -p $threads -ss $strand -1 ${filename}.fq $referencegenome $filename2"_SPRINT_output" bwa samtools
            fi
            rm $filename2.sam
        fi
    done

elif [ ! -z "$bam" ]; then

    for b in "${bam[@]}"; do
        if [ ! -z "$repeat" ]; then
            sprint_from_bam -cd $cluster_distance -csad1 $cluster_size_alu_ad1 -csad2 $cluster_size_alu_ad2 -csnar $cluster_size_nalurp -csnr $cluster_size_nrp -csrg $cluster_size_rg -rp $repeat "$b" $referencegenome $b"_SPRINT_output" samtools
        else
            sprint_from_bam -cd $cluster_distance -csrg $cluster_size_rg "$b" $referencegenome $b"_SPRINT_output" samtools
        fi   
    done
fi