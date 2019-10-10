#!/bin/bash

# set -x
# set -e

usage() {
      echo ""
      echo "Usage : sh $0 -g <reference_genome>  -i <Index_folder> -A <reference_annotation> -l lib_type {-1 <left_reads> -2 <right_reads> | -U <single_reads> | -s <sra_id>} -O <output_folder for Bam files> -p num_threads -5 <integer> -3 <integer> {-q phred_33 -Q phred_64} -m min_intron -M max_intron {-t stringtie -c cufflinks} -f <integer> -e <cuff_m>"
      echo ""

cat <<'EOF'
  
  ###### Command line options ##########

  -g <reference genome fasta file>

  -i <reference genomeindex folder>

  -A <reference genome annotation>

  -l Library type

  -1 <reads_1>
               # Make sure both the paired end reads are present
  -2 <reads_2>

  -U <single_reads> # Don not use Single Reads along with Paired end reads

  -O </path/to/bam output folder>

  -s SRA ID # One SRA ID or multiple SRA ID's in a file 

  -p Number of threads
  
  -5 5' trim 
 
  -3 3' trim

  -q phred33

  -Q phred64 

  -m Minimum intron length

  -M Maximum intron length

  -t StringTie (Report alignments tailored for transcript assemblers including StringTie)

  -c Cufflinks (Report alignments tailored specifically for Cufflinks)
  
  -f threshold # FPKM threshold to filter

  -e Cuffmerge

EOF
    exit 0
}

quality_33=0
quality_64=0
tra_as=0
tra_cuff=0
cuff_merge=0


while getopts ":hg:i:A:l:1:2:U:O:s:p:5:3:f:qQtcem:M:" opt; do
  case $opt in
    g)
    referencegenome=$OPTARG # Reference genome file
     ;;
    i)
    index_folder=$OPTARG # Input folder
     ;;
    A)
    referenceannotation=$OPTARG # Reference genome annotation
     ;;
    l)
     lib_type=$OPTARG # Library type
     ;;
    1)
    left_reads+=("$OPTARG") # Left reads
     ;;
    2)
    right_reads=("$OPTARG") # Right reads
     ;;
    U)
    single_reads+=("$OPTARG") # single end reads
     ;;
    O)
    bam_out=$OPTARG # Samoutput file
     ;;
    s)
    sra_id=$OPTARG # SRA ID or SRA ID's in a file
     ;;
    p)
    num_threads=$OPTARG # Number of threads
     ;;
    5)
    five_trim=$OPTARG # 5' trim
     ;;
    3)
    three_trim=$OPTARG # 3' trim
     ;;
    q)
    quality_33=$OPTARG # Phred 33 
     ;;
    Q)
    quality_64=$OPTARG # Phread 64
     ;;
    m)
    min_intl=$OPTARG # Minimum intron length
     ;;
    M)
    max_intl=$OPTARG # Maximum intron length
     ;;
    t)
    tra_as=$OPTARG # Report alignments tailored for transcript assemblers including StringTie
     ;;
    c)
    tra_cuff=$OPTARG # Report alignments tailored specifically for Cufflinks
     ;;
    f)
    threshold=$OPTARG # Coverage/base filter that you would like to apply to identified transcripts
     ;;
    e)
    cuff_merge=$OPTARG # Cuffmerge
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

# ############################################################################################################################################################################################################################
## Functions ###
# ############################################################################################################################################################################################################################

#Define parameters for coverage/base filtering
param5=5
param4=4
param3=3
param2=2
param1=1
param0=0

##################
### Single SRA ###
##################

coverge_cuffoff_SRA_single()
{
    if [ "$threshold" -eq "$param5" ]; then
       grep " transcript" "$bam_out"/$sra_id.gtf | grep -e 'cov "4.' -e 'cov "3.' -e 'cov "2.' -e 'cov "1.' -e 'cov "0.' | cut -f 9 | cut -d " " -f 4 | sort -u > listtoremove.txt
       grep -vFf listtoremove.txt "$bam_out"/$sra_id.gtf >"$bam_out"/$sra_id.gtf.filtered.gtf
       var=$(grep -vFf listtoremove.txt "$bam_out"/$sra_id.gtf.filtered.gtf | wc -l)
       if [ "$var" -eq 0 ]; then
          echo No transcripts have FPKM values exceeding your "$param5" cut-off have been found in $sra_id.gtf.filtered.gtf. Try lowering your cut-off 1>&2
       fi             
    elif [ "$threshold" -eq "$param4" ]; then
       grep " transcript" "$bam_out"/$sra_id.gtf | grep -e 'cov "3.' -e 'cov "2.' -e 'cov "1.' -e 'cov "0.' | cut -f 9 | cut -d " " -f 4 | sort -u > listtoremove.txt
       grep -vFf listtoremove.txt "$bam_out"/$sra_id.gtf >"$bam_out"/$sra_id.gtf.filtered.gtf
       var=$(grep -vFf listtoremove.txt "$bam_out"/$sra_id.gtf.filtered.gtf | wc -l)
       if [ "$var" -eq 0 ]; then
          echo No transcripts have FPKM values exceeding your "$param4" cut-off have been found in $sra_id.gtf.filtered.gtf. Try lowering your cut-off 1>&2
       fi             
    elif [ "$threshold" -eq "$param3" ]; then
       grep " transcript" "$bam_out"/$sra_id.gtf | grep -e 'cov "2.' -e 'cov "1.' -e 'cov "0.' | cut -f 9 | cut -d " " -f 4 | sort -u > listtoremove.txt
       grep -vFf listtoremove.txt "$bam_out"/$sra_id.gtf >"$bam_out"/$sra_id.gtf.filtered.gtf
       var=$(grep -vFf listtoremove.txt "$bam_out"/$sra_id.gtf.filtered.gtf | wc -l)
       if [ "$var" -eq 0 ]; then
          echo No transcripts have FPKM values exceeding your "$param3" cut-off have been found in $sra_id.gtf.filtered.gtf. Try lowering your cut-off 1>&2
       fi             
    elif [ "$threshold" -eq "$param2" ]; then
       grep " transcript" "$bam_out"/$sra_id.gtf | grep -e 'cov "1.' -e 'cov "0.' | cut -f 9 | cut -d " " -f 4 | sort -u > listtoremove.txt
       grep -vFf listtoremove.txt "$bam_out"/$sra_id.gtf >"$bam_out"/$sra_id.gtf.filtered.gtf
       var=$(grep -vFf listtoremove.txt "$bam_out"/$sra_id.gtf.filtered.gtf | wc -l)
       if [ "$var" -eq 0 ]; then
          echo No transcripts have FPKM values exceeding your "$param2" cut-off have been found in $sra_id.gtf.filtered.gtf. Try lowering your cut-off 1>&2
       fi             
    elif [ "$threshold" -eq "$param1" ]; then
       grep " transcript" "$bam_out"/$sra_id.gtf | grep -e 'cov "0.' | cut -f 9 | cut -d " " -f 4 | sort -u > listtoremove.txt
       grep -vFf listtoremove.txt "$bam_out"/$sra_id.gtf >"$bam_out"/$sra_id.gtf.filtered.gtf
       var=$(grep -vFf listtoremove.txt "$bam_out"/$sra_id.gtf.filtered.gtf | wc -l)
       if [ "$var" -eq 0 ]; then
          echo No transcripts have FPKM values exceeding your "$param1" cut-off have been found in $sra_id.gtf.filtered.gtf. Try lowering your cut-off 1>&2
       fi             
    elif [ "$threshold" -eq "$param0" ]; then
       grep " transcript" "$bam_out"/$sra_id.gtf | grep -e 'cov "0.000' | cut -f 9 | cut -d " " -f 4 | sort -u > listtoremove.txt
       grep -vFf listtoremove.txt "$bam_out"/$sra_id.gtf >"$bam_out"/$sra_id.gtf.filtered.gtf
    else
        echo "Invalid coverage parameter. Please select a whole number between 0-5"
        exit
    fi
}

stringtie_SRA_single()
{
    java -jar /usr/bin/picard.jar SortSam I=$sra_id.sam O="$bam_out"/$sra_id.sorted.bam SORT_ORDER=coordinate
    rm $sra_id.sam
    stringtie -G $referenceannotation "$bam_out"/$sra_id.sorted.bam -o "$bam_out"/$sra_id.gtf -p $num_threads
    coverge_cuffoff_SRA_single  
    var=$(grep -vFf listtoremove.txt "$bam_out"/$sra_id.gtf.filtered.gtf | wc -l)
    if [ "$var" -eq 0 ]; then
      rm listtoremove.txt "$bam_out"/$sra_id.gtf.filtered.gtf
    else
      cuffcompare "$bam_out"/$sra_id.gtf.filtered.gtf -r $referenceannotation -o $sra_id
      rm listtoremove.txt
      mv *.tracking *.loci *.combined.gtf *.stats "$bam_out"
   fi
}

cufflinks_SRA_single()
{
    java -jar /usr/bin/picard.jar SortSam I=$sra_id.sam O="$bam_out"/$sra_id.sorted.bam SORT_ORDER=coordinate
    rm $sra_id.sam
    cufflinks "$bam_out"/$sra_id.sorted.bam -p $num_threads -g $referenceannotation -o "$bam_out"
    mv "$bam_out"/transcripts.gtf "$bam_out"/$sra_id.gtf
    mv "$bam_out"/isoforms.fpkm_tracking "$bam_out"/$sra_id.isoforms.fpkm_tracking
    mv "$bam_out"/genes.fpkm_tracking "$bam_out"/$sra_id.genes.fpkm_tracking
    coverge_cuffoff_SRA_single           
    if [ "$var" -eq 0 ]; then
      rm listtoremove.txt "$bam_out"/$sra_id.gtf.filtered.gtf
    else
      cuffcompare "$bam_out"/$sra_id.gtf.filtered.gtf -r $referenceannotation -o $sra_id
      rm listtoremove.txt
      mv *.tracking *.loci *.combined.gtf *.stats "$bam_out"
   fi

}

##################
### Multi SRA ###
##################

coverge_cuffoff_SRA_multi()
{
    if [ "$threshold" -eq "$param5" ]; then
       grep " transcript" "$bam_out"/$f.gtf | grep -e 'cov "4.' -e 'cov "3.' -e 'cov "2.' -e 'cov "1.' -e 'cov "0.' | cut -f 9 | cut -d " " -f 4 | sort -u > listtoremove.txt
       grep -vFf listtoremove.txt "$bam_out"/$f.gtf >"$bam_out"/$f.gtf.filtered.gtf
       var=$(grep -vFf listtoremove.txt "$bam_out"/$f.gtf.filtered.gtf | wc -l)
       if [ "$var" -eq 0 ]; then
          echo No transcripts have FPKM values exceeding your "$param5" cut-off have been found in $f.gtf.filtered.gtf. Try lowering your cut-off 1>&2
       fi      
    elif [ "$threshold" -eq "$param4" ]; then
       grep " transcript" "$bam_out"/$f.gtf | grep -e 'cov "3.' -e 'cov "2.' -e 'cov "1.' -e 'cov "0.' | cut -f 9 | cut -d " " -f 4 | sort -u > listtoremove.txt
       grep -vFf listtoremove.txt "$bam_out"/$f.gtf >"$bam_out"/$f.gtf.filtered.gtf
       var=$(grep -vFf listtoremove.txt "$bam_out"/$f.gtf.filtered.gtf | wc -l)
       if [ "$var" -eq 0 ]; then
          echo No transcripts have FPKM values exceeding your "$param4" cut-off have been found in $f.gtf.filtered.gtf. Try lowering your cut-off 1>&2
       fi       
    elif [ "$threshold" -eq "$param3" ]; then
       grep " transcript" "$bam_out"/$f.gtf | grep -e 'cov "2.' -e 'cov "1.' -e 'cov "0.' | cut -f 9 | cut -d " " -f 4 | sort -u > listtoremove.txt
       grep -vFf listtoremove.txt "$bam_out"/$f.gtf >"$bam_out"/$f.gtf.filtered.gtf
       var=$(grep -vFf listtoremove.txt "$bam_out"/$f.gtf.filtered.gtf | wc -l)
       if [ "$var" -eq 0 ]; then
          echo No transcripts have FPKM values exceeding your "$param3" cut-off have been found in $f.gtf.filtered.gtf. Try lowering your cut-off 1>&2
       fi      
    elif [ "$threshold" -eq "$param2" ]; then
       grep " transcript" "$bam_out"/$f.gtf | grep -e 'cov "1.' -e 'cov "0.' | cut -f 9 | cut -d " " -f 4 | sort -u > listtoremove.txt
       grep -vFf listtoremove.txt "$bam_out"/$f.gtf >"$bam_out"/$f.gtf.filtered.gtf
       var=$(grep -vFf listtoremove.txt "$bam_out"/$f.gtf.filtered.gtf | wc -l)
       if [ "$var" -eq 0 ]; then
          echo No transcripts have FPKM values exceeding your "$param2" cut-off have been found in $f.gtf.filtered.gtf. Try lowering your cut-off 1>&2
       fi       
    elif [ "$threshold" -eq "$param1" ]; then
       grep " transcript" "$bam_out"/$f.gtf | grep -e 'cov "0.' | cut -f 9 | cut -d " " -f 4 | sort -u > listtoremove.txt
       grep -vFf listtoremove.txt "$bam_out"/$f.gtf >"$bam_out"/$f.gtf.filtered.gtf
       var=$(grep -vFf listtoremove.txt "$bam_out"/$f.gtf.filtered.gtf | wc -l)
       if [ "$var" -eq 0 ]; then
          echo No transcripts have FPKM values exceeding your "$param1" cut-off have been found in $f.gtf.filtered.gtf. Try lowering your cut-off 1>&2
       fi       
    elif [ "$threshold" -eq "$param0" ]; then
       grep " transcript" "$bam_out"/$f.gtf | grep -e 'cov "0.000' | cut -f 9 | cut -d " " -f 4 | sort -u > listtoremove.txt
       grep -vFf listtoremove.txt "$bam_out"/$f.gtf >"$bam_out"/$f.gtf.filtered.gtf
    else
        echo "Invalid coverage parameter. Please select a whole number between 0-5"
        exit
    fi
}

stringtie_SRA_multi()
{
    java -jar /usr/bin/picard.jar SortSam I=$f.sam O="$bam_out"/$f.sorted.bam SORT_ORDER=coordinate
    rm $f.sam
    stringtie -G $referenceannotation "$bam_out"/$f.sorted.bam -o "$bam_out"/$f.gtf -p $num_threads
    coverge_cuffoff_SRA_multi
    var=$(grep -vFf listtoremove.txt "$bam_out"/$f.gtf.filtered.gtf | wc -l)
    if [ "$var" -eq 0 ]; then
      rm listtoremove.txt "$bam_out"/$f.gtf.filtered.gtf
    else
      cuffcompare "$bam_out"/$f.gtf.filtered.gtf -r $referenceannotation -o $f
      rm listtoremove.txt
      mv *.tracking *.loci *.combined.gtf *.stats "$bam_out"
    fi
}

cufflinks_SRA_multi()
{
    java -jar /usr/bin/picard.jar SortSam I=$f.sam O="$bam_out"/$f.sorted.bam SORT_ORDER=coordinate
    rm $f.sam
    cufflinks "$bam_out"/$f.sorted.bam -p $num_threads -g $referenceannotation -o "$bam_out"
    mv "$bam_out"/transcripts.gtf "$bam_out"/$f.gtf
    mv "$bam_out"/isoforms.fpkm_tracking "$bam_out"/$f.isoforms.fpkm_tracking
    mv "$bam_out"/skipped.gtf  "$bam_out"/$f.skipped.gtf
    coverge_cuffoff_SRA_multi
    var=$(grep -vFf listtoremove.txt "$bam_out"/$f.gtf.filtered.gtf | wc -l)
    if [ "$var" -eq 0 ]; then
      rm listtoremove.txt "$bam_out"/$f.gtf.filtered.gtf
    else
      cuffcompare "$bam_out"/$f.gtf.filtered.gtf -r $referenceannotation -o $f
      rm listtoremove.txt
      mv *.tracking *.loci *.combined.gtf *.stats "$bam_out"
   fi

}

##################
### No SRA #######
##################

### Paired end #####

coverge_cuffoff_non_SRA()
{
    if [ "$threshold" -eq "$param5" ]; then
       grep " transcript" "$bam_out"/$filename3.gtf | grep -e 'cov "4.' -e 'cov "3.' -e 'cov "2.' -e 'cov "1.' -e 'cov "0.' | cut -f 9 | cut -d " " -f 4 | sort -u > listtoremove.txt
       grep -vFf listtoremove.txt "$bam_out"/$filename3.gtf >"$bam_out"/$filename3.gtf.filtered.gtf
       var=$(grep -vFf listtoremove.txt "$bam_out"/$filename3.gtf.filtered.gtf | wc -l)
       if [ "$var" -eq 0 ]; then
          echo No transcripts have FPKM values exceeding your "$param5" cut-off have been found in $filename3.gtf.filtered.gtf. Try lowering your cut-off 1>&2
       fi
    elif [ "$threshold" -eq "$param4" ]; then
       grep " transcript" "$bam_out"/$filename3.gtf | grep -e 'cov "3.' -e 'cov "2.' -e 'cov "1.' -e 'cov "0.' | cut -f 9 | cut -d " " -f 4 | sort -u > listtoremove.txt
       grep -vFf listtoremove.txt "$bam_out"/$filename3.gtf >"$bam_out"/$filename3.gtf.filtered.gtf
       var=$(grep -vFf listtoremove.txt "$bam_out"/$filename3.gtf.filtered.gtf | wc -l)
       if [ "$var" -eq 0 ]; then
          echo No transcripts have FPKM values exceeding your "$param4" cut-off have been found in $filename3.gtf.filtered.gtf. Try lowering your cut-off 1>&2
       fi
    elif [ "$threshold" -eq "$param3" ]; then
       grep " transcript" "$bam_out"/$filename3.gtf | grep -e 'cov "2.' -e 'cov "1.' -e 'cov "0.' | cut -f 9 | cut -d " " -f 4 | sort -u > listtoremove.txt
       grep -vFf listtoremove.txt "$bam_out"/$filename3.gtf >"$bam_out"/$filename3.gtf.filtered.gtf
       var=$(grep -vFf listtoremove.txt "$bam_out"/$filename3.gtf.filtered.gtf | wc -l)
       if [ "$var" -eq 0 ]; then
          echo No transcripts have FPKM values exceeding your "$param3" cut-off have been found in $filename3.gtf.filtered.gtf. Try lowering your cut-off 1>&2
       fi
    elif [ "$threshold" -eq "$param2" ]; then
       grep " transcript" "$bam_out"/$filename3.gtf | grep -e 'cov "1.' -e 'cov "0.' | cut -f 9 | cut -d " " -f 4 | sort -u > listtoremove.txt
       grep -vFf listtoremove.txt "$bam_out"/$filename3.gtf >"$bam_out"/$filename3.gtf.filtered.gtf
       var=$(grep -vFf listtoremove.txt "$bam_out"/$filename3.gtf.filtered.gtf | wc -l)
       if [ "$var" -eq 0 ]; then
          echo No transcripts have FPKM values exceeding your "$param2" cut-off have been found in $filename3.gtf.filtered.gtf. Try lowering your cut-off 1>&2
       fi
    elif [ "$threshold" -eq "$param1" ]; then
       grep " transcript" "$bam_out"/$filename3.gtf | grep -e 'cov "0.' | cut -f 9 | cut -d " " -f 4 | sort -u > listtoremove.txt
       grep -vFf listtoremove.txt "$bam_out"/$filename3.gtf >"$bam_out"/$filename3.gtf.filtered.gtf
       var=$(grep -vFf listtoremove.txt "$bam_out"/$filename3.gtf.filtered.gtf | wc -l)
       if [ "$var" -eq 0 ]; then
          echo No transcripts have FPKM values exceeding your "$param1" cut-off have been found in $filename3.gtf.filtered.gtf. Try lowering your cut-off 1>&2
       fi
    elif [ "$threshold" -eq "$param0" ]; then
       grep " transcript" "$bam_out"/$filename3.gtf | grep -e 'cov "0.000' | cut -f 9 | cut -d " " -f 4 | sort -u > listtoremove.txt
       grep -vFf listtoremove.txt "$bam_out"/$filename3.gtf >"$bam_out"/$filename3.gtf.filtered.gtf
    else
        echo "Invalid coverage parameter. Please select a whole number between 0-5"
        exit
    fi
}

stringtie_non_SRA() 
{
    java -jar /usr/bin/picard.jar SortSam I=$filename3.sam O="$bam_out"/$filename3.sorted.bam SORT_ORDER=coordinate
    rm $filename3.sam
    stringtie -G $referenceannotation "$bam_out"/$filename3.sorted.bam -o "$bam_out"/$filename3.gtf -p $num_threads
    coverge_cuffoff_non_SRA 
    var=$(grep -vFf listtoremove.txt "$bam_out"/$filename3.gtf.filtered.gtf | wc -l)
       if [ "$var" -eq 0 ]; then
          rm listtoremove.txt "$bam_out"/$filename3.gtf.filtered.gtf
        else
          cuffcompare "$bam_out"/$filename3.gtf.filtered.gtf -r $referenceannotation -o $filename3
          rm listtoremove.txt
          mv *.tracking *.loci *.combined.gtf *.stats "$bam_out"
       fi
}

cufflinks_non_SRA()
{
    java -jar /usr/bin/picard.jar SortSam I=$filename3.sam O="$bam_out"/$filename3.sorted.bam SORT_ORDER=coordinate
    rm $filename3.sam
    cufflinks "$bam_out"/$filename3.sorted.bam -p $num_threads -g $referenceannotation -o "$bam_out"
    mv "$bam_out"/skipped.gtf "$bam_out"/$filename3.skipped.gtf
    mv "$bam_out"/transcripts.gtf "$bam_out"/$filename3.gtf
    mv "$bam_out"/isoforms.fpkm_tracking "$bam_out"/$filename3.isoforms.fpkm_tracking
    mv "$bam_out"/genes.fpkm_tracking "$bam_out"/$filename3.genes.fpkm_tracking
    coverge_cuffoff_non_SRA
    var=$(grep -vFf listtoremove.txt "$bam_out"/$filename3.gtf.filtered.gtf | wc -l)
       if [ "$var" -eq 0 ]; then
          rm listtoremove.txt "$bam_out"/$filename3.gtf.filtered.gtf
        else
          cuffcompare "$bam_out"/$filename3.gtf.filtered.gtf -r $referenceannotation -o $filename3
          rm listtoremove.txt
          mv *.tracking *.loci *.combined.gtf *.stats "$bam_out"
       fi  
}

### single end ####

coverge_cuffoff_non_SRA_single()
{
    if [ "$threshold" -eq "$param5" ]; then
       grep " transcript" "$bam_out"/$filename.gtf | grep -e 'cov "4.' -e 'cov "3.' -e 'cov "2.' -e 'cov "1.' -e 'cov "0.' | cut -f 9 | cut -d " " -f 4 | sort -u > listtoremove.txt
       grep -vFf listtoremove.txt "$bam_out"/$filename.gtf >"$bam_out"/$filename.gtf.filtered.gtf
       var=$(grep -vFf listtoremove.txt "$bam_out"/$filename.gtf.filtered.gtf | wc -l)
       if [ "$var" -eq 0 ]; then
          echo No transcripts have FPKM values exceeding your "$param5" cut-off have been found in $filename.gtf.filtered.gtf. Try lowering your cut-off 1>&2
       fi       
    elif [ "$threshold" -eq "$param4" ]; then
       grep " transcript" "$bam_out"/$filename.gtf | grep -e 'cov "3.' -e 'cov "2.' -e 'cov "1.' -e 'cov "0.' | cut -f 9 | cut -d " " -f 4 | sort -u > listtoremove.txt
       grep -vFf listtoremove.txt "$bam_out"/$filename.gtf >"$bam_out"/$filename.gtf.filtered.gtf
       var=$(grep -vFf listtoremove.txt "$bam_out"/$filename.gtf.filtered.gtf | wc -l)
       if [ "$var" -eq 0 ]; then
          echo No transcripts have FPKM values exceeding your "$param4" cut-off have been found in $filename.gtf.filtered.gtf. Try lowering your cut-off 1>&2
       fi
    elif [ "$threshold" -eq "$param3" ]; then
       grep " transcript" "$bam_out"/$filename.gtf | grep -e 'cov "2.' -e 'cov "1.' -e 'cov "0.' | cut -f 9 | cut -d " " -f 4 | sort -u > listtoremove.txt
       grep -vFf listtoremove.txt "$bam_out"/$filename.gtf >"$bam_out"/$filename.gtf.filtered.gtf
       var=$(grep -vFf listtoremove.txt "$bam_out"/$filename.gtf.filtered.gtf | wc -l)
       if [ "$var" -eq 0 ]; then
          echo No transcripts have FPKM values exceeding your "$param3" cut-off have been found in $filename.gtf.filtered.gtf. Try lowering your cut-off 1>&2
       fi
    elif [ "$threshold" -eq "$param2" ]; then
       grep " transcript" "$bam_out"/$filename.gtf | grep -e 'cov "1.' -e 'cov "0.' | cut -f 9 | cut -d " " -f 4 | sort -u > listtoremove.txt
       grep -vFf listtoremove.txt "$bam_out"/$filename.gtf >"$bam_out"/$filename.gtf.filtered.gtf
       var=$(grep -vFf listtoremove.txt "$bam_out"/$filename.gtf.filtered.gtf | wc -l)
       if [ "$var" -eq 0 ]; then
          echo No transcripts have FPKM values exceeding your "$param2" cut-off have been found in $filename.gtf.filtered.gtf. Try lowering your cut-off 1>&2
       fi
    elif [ "$threshold" -eq "$param1" ]; then
       grep " transcript" "$bam_out"/$filename.gtf | grep -e 'cov "0.' | cut -f 9 | cut -d " " -f 4 | sort -u > listtoremove.txt
       grep -vFf listtoremove.txt "$bam_out"/$filename.gtf >"$bam_out"/$filename.gtf.filtered.gtf
       var=$(grep -vFf listtoremove.txt "$bam_out"/$filename.gtf.filtered.gtf | wc -l)
       if [ "$var" -eq 0 ]; then
          echo No transcripts have FPKM values exceeding your "$param1" cut-off have been found in $filename.gtf.filtered.gtf. Try lowering your cut-off 1>&2
       fi
    elif [ "$threshold" -eq "$param0" ]; then
       grep " transcript" "$bam_out"/$filename.gtf | grep -e 'cov "0.000' | cut -f 9 | cut -d " " -f 4 | sort -u > listtoremove.txt
       grep -vFf listtoremove.txt "$bam_out"/$filename.gtf >"$bam_out"/$filename.gtf.filtered.gtf
    else
        echo "Invalid coverage parameter. Please select a whole number between 0-5"
        exit
    fi
}
 

stringtie_non_SRA_single() 
{
    java -jar /usr/bin/picard.jar SortSam I=$filename.sam O="$bam_out"/$filename.sorted.bam SORT_ORDER=coordinate
    rm $filename.sam
    stringtie -G $referenceannotation "$bam_out"/$filename.sorted.bam -o "$bam_out"/$filename.gtf -p $num_threads
    coverge_cuffoff_non_SRA_single
    var=$(grep -vFf listtoremove.txt "$bam_out"/$filename.gtf.filtered.gtf | wc -l)
       if [ "$var" -eq 0 ]; then
          rm listtoremove.txt "$bam_out"/$filename.gtf.filtered.gtf
        else
          cuffcompare "$bam_out"/$filename.gtf.filtered.gtf -r $referenceannotation -o $filename
          rm listtoremove.txt
          mv *.tracking *.loci *.combined.gtf *.stats "$bam_out"
       fi
}

cufflinks_non_SRA_single()
{
    java -jar /usr/bin/picard.jar SortSam I=$filename.sam O="$bam_out"/$filename.sorted.bam SORT_ORDER=coordinate
    rm $filename.sam
    cufflinks "$bam_out"/$filename.sorted.bam -p $num_threads -g $referenceannotation -o "$bam_out"
    mv "$bam_out"/skipped.gtf "$bam_out"/$filename.skipped.gtf
    mv "$bam_out"/transcripts.gtf "$bam_out"/$filename.gtf
    mv "$bam_out"/isoforms.fpkm_tracking "$bam_out"/$filename.isoforms.fpkm_tracking
    mv "$bam_out"/genes.fpkm_tracking "$bam_out"/$filename.genes.fpkm_tracking
    coverge_cuffoff_non_SRA_single        
    var=$(grep -vFf listtoremove.txt "$bam_out"/$filename.gtf.filtered.gtf | wc -l)
       if [ "$var" -eq 0 ]; then
          rm listtoremove.txt "$bam_out"/$filename.gtf.filtered.gtf
        else
          cuffcompare "$bam_out"/$filename.gtf.filtered.gtf -r $referenceannotation -o $filename
          rm listtoremove.txt
          mv *.tracking *.loci *.combined.gtf *.stats "$bam_out"
       fi
}


###############
## cuffmerge ##
###############

cuff_merge_fun()
{
    files="$bam_out"/*filtered.gtf
    if [ ${#files[@]} -gt 0 ]; then
      ls "$bam_out"/*filtered.gtf | tr "\t" "\n" >> "$bam_out"/gtf_file.txt
      numb=$(cat "$bam_out"/gtf_file.txt | wc -l)
      if [ "$numb" -gt 1 ]; then  
        cuffmerge -o "$bam_out"/merged_out -g $referenceannotation "$bam_out"/gtf_file.txt -p $num_threads
        "$bam_out"/gtf_file.txt
      else
        echo "Cuffmerge needs more than 1 filtered.gtf file" 1>&2
      fi
    else
      echo "No filtered.gtf files are found for cuffmerge" 1>&2
    fi  
}

# ############################################################################################################################################################################################################################
# # Reference genome/Index
# ############################################################################################################################################################################################################################
# Index folder

if [ ! -z "$index_folder" ]; then
  for i in $index_folder/*; do
      cp $i .
      fbname=$(basename "$i" .ht2 | cut -d. -f1)
  done

elif [ ! -z "$referencegenome" ] && [ -z "$index_folder" ]; then
  hisat2-build $referencegenome ref_genome -p $num_threads
  fbname=$(basename "ref_genome" .ht2 | cut -d. -f1)
fi

# ############################################################################################################################################################################################################################
# # Paired end reads
# ############################################################################################################################################################################################################################

#### Stringtie #####

# Phred 33

if [ ! -z "$left_reads" ] && [ ! -z "$right_reads" ] && [ "$quality_33" != 0 ] && [ "$tra_as" != 0 ] && [ "$tra_cuff" == 0 ]; then
    mkdir "$bam_out"
    numb=$(ls "${left_reads[@]}" | wc -l)
    for f in "${left_reads[@]}"; do
      extension=$(echo "$f" | sed -r 's/.*(fq|fq.gz|fastq|fastq.gz)$/\1/')
       
      if [[ "$extension" =~ "fq.gz" ]]; then
        filename=$(basename "$f" ".fq.gz")
        filename2=${filename/_R1/_R2}
        filename3=$(echo $filename | sed 's/_R1//')
        
        if [ $numb -eq 1 ] && [ "$cuff_merge" == 0 ]; then
            hisat2 -x $fbname --rna-strandness $lib_type -1 ${filename}.fq.gz -2 ${filename2}.fq.gz -S $filename3.sam -p $num_threads -5 $five_trim -3 $three_trim --min-intronlen $min_intl --max-intronlen $max_intl --dta --phred33
            stringtie_non_SRA
            rm $fbname*
            exit
        elif [ $numb -eq 1 ] && [ "$cuff_merge" != 0 ]; then
            hisat2 -x $fbname --rna-strandness $lib_type -1 ${filename}.fq.gz -2 ${filename2}.fq.gz -S $filename3.sam -p $num_threads -5 $five_trim -3 $three_trim --min-intronlen $min_intl --max-intronlen $max_intl --dta --phred33
            stringtie_non_SRA
            rm $fbname*
            echo "cuffmerge only works with more than one file"
            exit
        elif [ $numb -gt 1 ] && [ "$cuff_merge" == 0 ]; then
            hisat2 -x $fbname --rna-strandness $lib_type -1 ${filename}.fq.gz -2 ${filename2}.fq.gz -S $filename3.sam -p $num_threads -5 $five_trim -3 $three_trim --min-intronlen $min_intl --max-intronlen $max_intl --dta --phred33
            stringtie_non_SRA   
        elif [ $numb -gt 1 ] && [ "$cuff_merge" != 0 ]; then
            hisat2 -x $fbname --rna-strandness $lib_type -1 ${filename}.fq.gz -2 ${filename2}.fq.gz -S $filename3.sam -p $num_threads -5 $five_trim -3 $three_trim --min-intronlen $min_intl --max-intronlen $max_intl --dta --phred33
            stringtie_non_SRA
        fi

      elif [[ "$extension" =~ "fastq.gz" ]]; then
        filename=$(basename "$f" ".fastq.gz")
        filename2=${filename/_R1/_R2}
	      filename3=$(echo $filename | sed 's/_R1//')

        if [ $numb -eq 1 ] && [ "$cuff_merge" == 0 ]; then
            hisat2 -x $fbname --rna-strandness $lib_type -1 ${filename}.fastq.gz -2 ${filename2}.fastq.gz -S $filename3.sam -p $num_threads -5 $five_trim -3 $three_trim --min-intronlen $min_intl --max-intronlen $max_intl --dta --phred33
            stringtie_non_SRA
            rm $fbname*
            exit
        elif [ $numb -eq 1 ] && [ "$cuff_merge" != 0 ]; then
            hisat2 -x $fbname --rna-strandness $lib_type -1 ${filename}.fastq.gz -2 ${filename2}.fastq.gz -S $filename3.sam -p $num_threads -5 $five_trim -3 $three_trim --min-intronlen $min_intl --max-intronlen $max_intl --dta --phred33
            stringtie_non_SRA
            rm $fbname*
            echo "cuffmerge only works with more than one file"
            exit    
        elif [ $numb -gt 1 ] && [ "$cuff_merge" == 0 ]; then
            hisat2 -x $fbname --rna-strandness $lib_type -1 ${filename}.fastq.gz -2 ${filename2}.fastq.gz -S $filename3.sam -p $num_threads -5 $five_trim -3 $three_trim --min-intronlen $min_intl --max-intronlen $max_intl --dta --phred33
            stringtie_non_SRA
        elif [ $numb -gt 1 ] && [ "$cuff_merge" != 0 ]; then
            hisat2 -x $fbname --rna-strandness $lib_type -1 ${filename}.fastq.gz -2 ${filename2}.fastq.gz -S $filename3.sam -p $num_threads -5 $five_trim -3 $three_trim --min-intronlen $min_intl --max-intronlen $max_intl --dta --phred33
            stringtie_non_SRA
        fi

      elif [[ "$extension" =~ "fq" ]]; then
        filename=$(basename "$f" ".fq")
        filename2=${filename/_R1/_R2}
	      filename3=$(echo $filename | sed 's/_R1//')

        if [ $numb -eq 1 ] && [ "$cuff_merge" == 0 ]; then
            hisat2 -x $fbname --rna-strandness $lib_type -1 ${filename}.fq -2 ${filename2}.fq -S $filename3.sam -p $num_threads -5 $five_trim -3 $three_trim --min-intronlen $min_intl --max-intronlen $max_intl --dta --phred33
            stringtie_non_SRA
            rm $fbname*
            exit
        elif [ $numb -eq 1 ] && [ "$cuff_merge" != 0 ]; then
            hisat2 -x $fbname --rna-strandness $lib_type -1 ${filename}.fq -2 ${filename2}.fq -S $filename3.sam -p $num_threads -5 $five_trim -3 $three_trim --min-intronlen $min_intl --max-intronlen $max_intl --dta --phred33
            stringtie_non_SRA
            rm $fbname*
            echo "cuffmerge only works with more than one file"
            exit
        elif [ $numb -gt 1 ] && [ "$cuff_merge" == 0 ]; then
            hisat2 -x $fbname --rna-strandness $lib_type -1 ${filename}.fq -2 ${filename2}.fq -S $filename3.sam -p $num_threads -5 $five_trim -3 $three_trim --min-intronlen $min_intl --max-intronlen $max_intl --dta --phred33
            stringtie_non_SRA   
        elif [ $numb -gt 1 ] && [ "$cuff_merge" != 0 ]; then
            hisat2 -x $fbname --rna-strandness $lib_type -1 ${filename}.fq -2 ${filename2}.fq -S $filename3.sam -p $num_threads -5 $five_trim -3 $three_trim --min-intronlen $min_intl --max-intronlen $max_intl --dta --phred33
            stringtie_non_SRA
        fi

      elif [[ "$extension" =~ "fastq" ]]; then
        filename=$(basename "$f" ".fastq")
        filename2=${filename/_R1/_R2}
	      filename3=$(echo $filename | sed 's/_R1//')

        if [ $numb -eq 1 ] && [ "$cuff_merge" == 0 ]; then
            hisat2 -x $fbname --rna-strandness $lib_type -1 ${filename}.fastq -2 ${filename2}.fastq -S $filename3.sam -p $num_threads -5 $five_trim -3 $three_trim --min-intronlen $min_intl --max-intronlen $max_intl --dta --phred33
            stringtie_non_SRA
            rm $fbname*
            exit
        elif [ $numb -eq 1 ] && [ "$cuff_merge" != 0 ]; then
            hisat2 -x $fbname --rna-strandness $lib_type -1 ${filename}.fastq -2 ${filename2}.fastq -S $filename3.sam -p $num_threads -5 $five_trim -3 $three_trim --min-intronlen $min_intl --max-intronlen $max_intl --dta --phred33
            stringtie_non_SRA
            rm $fbname*
            echo "cuffmerge only works with more than one file"
            exit
        elif [ $numb -gt 1 ] && [ "$cuff_merge" == 0 ]; then
            hisat2 -x $fbname --rna-strandness $lib_type -1 ${filename}.fastq -2 ${filename2}.fastq -S $filename3.sam -p $num_threads -5 $five_trim -3 $three_trim --min-intronlen $min_intl --max-intronlen $max_intl --dta --phred33
            stringtie_non_SRA               
        elif [ $numb -gt 1 ] && [ "$cuff_merge" != 0 ]; then
            hisat2 -x $fbname --rna-strandness $lib_type -1 ${filename}.fastq -2 ${filename2}.fastq -S $filename3.sam -p $num_threads -5 $five_trim -3 $three_trim --min-intronlen $min_intl --max-intronlen $max_intl --dta --phred33
            stringtie_non_SRA
        fi   

      elif [ "$extension" != "fastq" ] || [ "$extension" != "fq" ] || [ "$extension" != "fastq.gz" ] || [ "$extension" != "fq.gz" ]; then
        rm $fbname*
        echo "The extension" "$extension" "is not supported. Only .fq, .fq.gz, .fastq, .fastq.gz are only supported" 1>&2        
        exit 64
      fi 
    done
    
    if [ $numb -gt 1 ] && [ "$cuff_merge" == 0 ]; then
      exit
    elif [ $numb -gt 1 ] && [ "$cuff_merge" != 0 ]; then
      cuff_merge_fun
    fi

# Phred 64

elif [ ! -z "$left_reads" ] && [ ! -z "$right_reads" ] && [ "$quality_64" != 0 ] && [ "$tra_as" != 0 ] && [ "$tra_cuff" == 0 ]; then
    mkdir "$bam_out"
    numb=$(ls "${left_reads[@]}" | wc -l)
    for f in "${left_reads[@]}"; do
      extension=$(echo "$f" | sed -r 's/.*(fq|fq.gz|fastq|fastq.gz)$/\1/')
       
      if [[ "$extension" =~ "fq.gz" ]]; then
        filename=$(basename "$f" ".fq.gz")
        filename2=${filename/_R1/_R2}
        filename3=$(echo $filename | sed 's/_R1//')

        if [ $numb -eq 1 ] && [ "$cuff_merge" == 0 ]; then
            hisat2 -x $fbname --rna-strandness $lib_type -1 ${filename}.fq.gz -2 ${filename2}.fq.gz -S $filename3.sam -p $num_threads -5 $five_trim -3 $three_trim --min-intronlen $min_intl --max-intronlen $max_intl --dta --phred64
            stringtie_non_SRA
            rm $fbname*
            exit
        elif [ $numb -eq 1 ] && [ "$cuff_merge" != 0 ]; then
            hisat2 -x $fbname --rna-strandness $lib_type -1 ${filename}.fq.gz -2 ${filename2}.fq.gz -S $filename3.sam -p $num_threads -5 $five_trim -3 $three_trim --min-intronlen $min_intl --max-intronlen $max_intl --dta --phred64
            stringtie_non_SRA
            rm $fbname*
            echo "cuffmerge only works with more than one file"
            exit    
        elif [ $numb -gt 1 ] && [ "$cuff_merge" == 0 ]; then
            hisat2 -x $fbname --rna-strandness $lib_type -1 ${filename}.fq.gz -2 ${filename2}.fq.gz -S $filename3.sam -p $num_threads -5 $five_trim -3 $three_trim --min-intronlen $min_intl --max-intronlen $max_intl --dta --phred64
            stringtie_non_SRA
        elif [ $numb -gt 1 ] && [ "$cuff_merge" != 0 ]; then
            hisat2 -x $fbname --rna-strandness $lib_type -1 ${filename}.fq.gz -2 ${filename2}.fq.gz -S $filename3.sam -p $num_threads -5 $five_trim -3 $three_trim --min-intronlen $min_intl --max-intronlen $max_intl --dta --phred64
            stringtie_non_SRA
        fi

      elif [[ "$extension" =~ "fastq.gz" ]]; then
        filename=$(basename "$f" ".fastq.gz")
        filename2=${filename/_R1/_R2}
	      filename3=$(echo $filename | sed 's/_R1//')

        if [ $numb -eq 1 ] && [ "$cuff_merge" == 0 ]; then
            hisat2 -x $fbname --rna-strandness $lib_type -1 ${filename}.fastq.gz -2 ${filename2}.fastq.gz -S $filename3.sam -p $num_threads -5 $five_trim -3 $three_trim --min-intronlen $min_intl --max-intronlen $max_intl --dta --phred64
            stringtie_non_SRA
            rm $fbname*
            exit
        elif [ $numb -eq 1 ] && [ "$cuff_merge" != 0 ]; then
            hisat2 -x $fbname --rna-strandness $lib_type -1 ${filename}.fastq.gz -2 ${filename2}.fastq.gz -S $filename3.sam -p $num_threads -5 $five_trim -3 $three_trim --min-intronlen $min_intl --max-intronlen $max_intl --dta --phred64
            stringtie_non_SRA
            rm $fbname*
            echo "cuffmerge only works with more than one file"
            exit    
        elif [ $numb -gt 1 ] && [ "$cuff_merge" == 0 ]; then
            hisat2 -x $fbname --rna-strandness $lib_type -1 ${filename}.fastq.gz -2 ${filename2}.fastq.gz -S $filename3.sam -p $num_threads -5 $five_trim -3 $three_trim --min-intronlen $min_intl --max-intronlen $max_intl --dta --phred64
            stringtie_non_SRA
        elif [ $numb -gt 1 ] && [ "$cuff_merge" != 0 ]; then
            hisat2 -x $fbname --rna-strandness $lib_type -1 ${filename}.fastq.gz -2 ${filename2}.fastq.gz -S $filename3.sam -p $num_threads -5 $five_trim -3 $three_trim --min-intronlen $min_intl --max-intronlen $max_intl --dta --phred64
            stringtie_non_SRA
        fi

      elif [[ "$extension" =~ "fq" ]]; then
        filename=$(basename "$f" ".fq")
        filename2=${filename/_R1/_R2}
	      filename3=$(echo $filename | sed 's/_R1//')

        if [ $numb -eq 1 ] && [ "$cuff_merge" == 0 ]; then
            hisat2 -x $fbname --rna-strandness $lib_type -1 ${filename}.fq -2 ${filename2}.fq -S $filename3.sam -p $num_threads -5 $five_trim -3 $three_trim --min-intronlen $min_intl --max-intronlen $max_intl --dta --phred64
            stringtie_non_SRA
            rm $fbname*
            exit
        elif [ $numb -eq 1 ] && [ "$cuff_merge" != 0 ]; then
            hisat2 -x $fbname --rna-strandness $lib_type -1 ${filename}.fq -2 ${filename2}.fq -S $filename3.sam -p $num_threads -5 $five_trim -3 $three_trim --min-intronlen $min_intl --max-intronlen $max_intl --dta --phred64
            stringtie_non_SRA
            rm $fbname*
            echo "cuffmerge only works with more than one file"
            exit    
        elif [ $numb -gt 1 ] && [ "$cuff_merge" == 0 ]; then
            hisat2 -x $fbname --rna-strandness $lib_type -1 ${filename}.fq -2 ${filename2}.fq -S $filename3.sam -p $num_threads -5 $five_trim -3 $three_trim --min-intronlen $min_intl --max-intronlen $max_intl --dta --phred64
            stringtie_non_SRA
        elif [ $numb -gt 1 ] && [ "$cuff_merge" != 0 ]; then
            hisat2 -x $fbname --rna-strandness $lib_type -1 ${filename}.fq -2 ${filename2}.fq -S $filename3.sam -p $num_threads -5 $five_trim -3 $three_trim --min-intronlen $min_intl --max-intronlen $max_intl --dta --phred64
            stringtie_non_SRA
        fi

      elif [[ "$extension" =~ "fastq" ]]; then
        filename=$(basename "$f" ".fastq")
        filename2=${filename/_R1/_R2}
	      filename3=$(echo $filename | sed 's/_R1//')

        if [ $numb -eq 1 ] && [ "$cuff_merge" == 0 ]; then
            hisat2 -x $fbname --rna-strandness $lib_type -1 ${filename}.fastq -2 ${filename2}.fastq -S $filename3.sam -p $num_threads -5 $five_trim -3 $three_trim --min-intronlen $min_intl --max-intronlen $max_intl --dta --phred64
            stringtie_non_SRA
            rm $fbname*
            exit
        elif [ $numb -eq 1 ] && [ "$cuff_merge" != 0 ]; then
            hisat2 -x $fbname --rna-strandness $lib_type -1 ${filename}.fastq -2 ${filename2}.fastq -S $filename3.sam -p $num_threads -5 $five_trim -3 $three_trim --min-intronlen $min_intl --max-intronlen $max_intl --dta --phred64
            stringtie_non_SRA
            rm $fbname*
            echo "cuffmerge only works with more than one file"
            exit    
        elif [ $numb -gt 1 ] && [ "$cuff_merge" == 0 ]; then
            hisat2 -x $fbname --rna-strandness $lib_type -1 ${filename}.fastq -2 ${filename2}.fastq -S $filename3.sam -p $num_threads -5 $five_trim -3 $three_trim --min-intronlen $min_intl --max-intronlen $max_intl --dta --phred64
            stringtie_non_SRA
        elif [ $numb -gt 1 ] && [ "$cuff_merge" != 0 ]; then
            hisat2 -x $fbname --rna-strandness $lib_type -1 ${filename}.fastq -2 ${filename2}.fastq -S $filename3.sam -p $num_threads -5 $five_trim -3 $three_trim --min-intronlen $min_intl --max-intronlen $max_intl --dta --phred64
            stringtie_non_SRA
        fi   

      elif [ "$extension" != "fastq" ] || [ "$extension" != "fq" ] || [ "$extension" != "fastq.gz" ] || [ "$extension" != "fq.gz" ]; then
        rm $fbname*
        echo "The extension" "$extension" "is not supported. Only .fq, .fq.gz, .fastq, .fastq.gz are only supported" 1>&2        
        exit 64
      fi 
    done
    
    if [ $numb -gt 1 ] && [ "$cuff_merge" == 0 ]; then
      exit
    elif [ $numb -gt 1 ] && [ "$cuff_merge" != 0 ]; then
      cuff_merge_fun
    fi


#### Cufflinks ######

# Phread-33

elif [ ! -z "$left_reads" ] && [ ! -z "$right_reads" ] && [ "$quality_33" != 0 ] && [ "$tra_as" == 0 ] && [ "$tra_cuff" != 0 ]; then
    mkdir "$bam_out"
    numb=$(ls "${left_reads[@]}" | wc -l)
    for f in "${left_reads[@]}"; do
      extension=$(echo "$f" | sed -r 's/.*(fq|fq.gz|fastq|fastq.gz)$/\1/')
       
      if [[ "$extension" =~ "fq.gz" ]]; then
        filename=$(basename "$f" ".fq.gz")
        filename2=${filename/_R1/_R2}
        filename3=$(echo $filename | sed 's/_R1//')

        if [ $numb -eq 1 ] && [ "$cuff_merge" == 0 ]; then
            hisat2 -x $fbname --rna-strandness $lib_type -1 ${filename}.fq.gz -2 ${filename2}.fq.gz -S $filename3.sam -p $num_threads -5 $five_trim -3 $three_trim --min-intronlen $min_intl --max-intronlen $max_intl --dta --phred33
            cufflinks_non_SRA
            rm $fbname*
            exit
        elif [ $numb -eq 1 ] && [ "$cuff_merge" != 0 ]; then
            hisat2 -x $fbname --rna-strandness $lib_type -1 ${filename}.fq.gz -2 ${filename2}.fq.gz -S $filename3.sam -p $num_threads -5 $five_trim -3 $three_trim --min-intronlen $min_intl --max-intronlen $max_intl --dta --phred33
            cufflinks_non_SRA
            rm $fbname*
            echo "cuffmerge only works with more than one file"
            exit    
        elif [ $numb -gt 1 ] && [ "$cuff_merge" == 0 ]; then
            hisat2 -x $fbname --rna-strandness $lib_type -1 ${filename}.fq.gz -2 ${filename2}.fq.gz -S $filename3.sam -p $num_threads -5 $five_trim -3 $three_trim --min-intronlen $min_intl --max-intronlen $max_intl --dta --phred33
            cufflinks_non_SRA
        elif [ $numb -gt 1 ] && [ "$cuff_merge" != 0 ]; then
            hisat2 -x $fbname --rna-strandness $lib_type -1 ${filename}.fq.gz -2 ${filename2}.fq.gz -S $filename3.sam -p $num_threads -5 $five_trim -3 $three_trim --min-intronlen $min_intl --max-intronlen $max_intl --dta --phred33
            cufflinks_non_SRA
        fi

      elif [[ "$extension" =~ "fastq.gz" ]]; then
        filename=$(basename "$f" ".fastq.gz")
        filename2=${filename/_R1/_R2}
	      filename3=$(echo $filename | sed 's/_R1//')

        if [ $numb -eq 1 ] && [ "$cuff_merge" == 0 ]; then
            hisat2 -x $fbname --rna-strandness $lib_type -1 ${filename}.fastq.gz -2 ${filename2}.fastq.gz -S $filename3.sam -p $num_threads -5 $five_trim -3 $three_trim --min-intronlen $min_intl --max-intronlen $max_intl --dta --phred33
            cufflinks_non_SRA
            rm $fbname*
            exit
        elif [ $numb -eq 1 ] && [ "$cuff_merge" != 0 ]; then
            hisat2 -x $fbname --rna-strandness $lib_type -1 ${filename}.fastq.gz -2 ${filename2}.fastq.gz -S $filename3.sam -p $num_threads -5 $five_trim -3 $three_trim --min-intronlen $min_intl --max-intronlen $max_intl --dta --phred33
            cufflinks_non_SRA
            rm $fbname*
            echo "cuffmerge only works with more than one file"
            exit    
        elif [ $numb -gt 1 ] && [ "$cuff_merge" == 0 ]; then
            hisat2 -x $fbname --rna-strandness $lib_type -1 ${filename}.fastq.gz -2 ${filename2}.fastq.gz -S $filename3.sam -p $num_threads -5 $five_trim -3 $three_trim --min-intronlen $min_intl --max-intronlen $max_intl --dta --phred33
            cufflinks_non_SRA
        elif [ $numb -gt 1 ] && [ "$cuff_merge" != 0 ]; then
            hisat2 -x $fbname --rna-strandness $lib_type -1 ${filename}.fastq.gz -2 ${filename2}.fastq.gz -S $filename3.sam -p $num_threads -5 $five_trim -3 $three_trim --min-intronlen $min_intl --max-intronlen $max_intl --dta --phred33
            cufflinks_non_SRA
        fi

      elif [[ "$extension" =~ "fq" ]]; then
        filename=$(basename "$f" ".fq")
        filename2=${filename/_R1/_R2}
	      filename3=$(echo $filename | sed 's/_R1//')

        if [ $numb -eq 1 ] && [ "$cuff_merge" == 0 ]; then
            hisat2 -x $fbname --rna-strandness $lib_type -1 ${filename}.fq -2 ${filename2}.fq -S $filename3.sam -p $num_threads -5 $five_trim -3 $three_trim --min-intronlen $min_intl --max-intronlen $max_intl --dta --phred33
            cufflinks_non_SRA
            rm $fbname*
            exit
        elif [ $numb -eq 1 ] && [ "$cuff_merge" != 0 ]; then
            hisat2 -x $fbname --rna-strandness $lib_type -1 ${filename}.fq -2 ${filename2}.fq -S $filename3.sam -p $num_threads -5 $five_trim -3 $three_trim --min-intronlen $min_intl --max-intronlen $max_intl --dta --phred33
            cufflinks_non_SRA
            rm $fbname*
            echo "cuffmerge only works with more than one file"
            exit    
        elif [ $numb -gt 1 ] && [ "$cuff_merge" == 0 ]; then
            hisat2 -x $fbname --rna-strandness $lib_type -1 ${filename}.fq -2 ${filename2}.fq -S $filename3.sam -p $num_threads -5 $five_trim -3 $three_trim --min-intronlen $min_intl --max-intronlen $max_intl --dta --phred33
            cufflinks_non_SRA
        elif [ $numb -gt 1 ] && [ "$cuff_merge" != 0 ]; then
            hisat2 -x $fbname --rna-strandness $lib_type -1 ${filename}.fq -2 ${filename2}.fq -S $filename3.sam -p $num_threads -5 $five_trim -3 $three_trim --min-intronlen $min_intl --max-intronlen $max_intl --dta --phred33
            cufflinks_non_SRA
        fi

      elif [[ "$extension" =~ "fastq" ]]; then
        filename=$(basename "$f" ".fastq")
        filename2=${filename/_R1/_R2}
	      filename3=$(echo $filename | sed 's/_R1//')

        if [ $numb -eq 1 ] && [ "$cuff_merge" == 0 ]; then
            hisat2 -x $fbname --rna-strandness $lib_type -1 ${filename}.fastq -2 ${filename2}.fastq -S $filename3.sam -p $num_threads -5 $five_trim -3 $three_trim --min-intronlen $min_intl --max-intronlen $max_intl --dta --phred33
            cufflinks_non_SRA
            rm $fbname*
            exit
        elif [ $numb -eq 1 ] && [ "$cuff_merge" != 0 ]; then
            hisat2 -x $fbname --rna-strandness $lib_type -1 ${filename}.fastq -2 ${filename2}.fastq -S $filename3.sam -p $num_threads -5 $five_trim -3 $three_trim --min-intronlen $min_intl --max-intronlen $max_intl --dta --phred33
            cufflinks_non_SRA
            rm $fbname*
            echo "cuffmerge only works with more than one file"
            exit    
        elif [ $numb -gt 1 ] && [ "$cuff_merge" == 0 ]; then
            hisat2 -x $fbname --rna-strandness $lib_type -1 ${filename}.fastq -2 ${filename2}.fastq -S $filename3.sam -p $num_threads -5 $five_trim -3 $three_trim --min-intronlen $min_intl --max-intronlen $max_intl --dta --phred33
            cufflinks_non_SRA
        elif [ $numb -gt 1 ] && [ "$cuff_merge" != 0 ]; then
            hisat2 -x $fbname --rna-strandness $lib_type -1 ${filename}.fastq -2 ${filename2}.fastq -S $filename3.sam -p $num_threads -5 $five_trim -3 $three_trim --min-intronlen $min_intl --max-intronlen $max_intl --dta --phred33
            cufflinks_non_SRA
        fi   

      elif [ "$extension" != "fastq" ] || [ "$extension" != "fq" ] || [ "$extension" != "fastq.gz" ] || [ "$extension" != "fq.gz" ]; then
        rm $fbname*
        echo "The extension" "$extension" "is not supported. Only .fq, .fq.gz, .fastq, .fastq.gz are only supported" 1>&2        
        exit 64
      fi 
    done
    
    if [ $numb -gt 1 ] && [ "$cuff_merge" == 0 ]; then
      exit
    elif [ $numb -gt 1 ] && [ "$cuff_merge" != 0 ]; then
      cuff_merge_fun
    fi


# Phred 64

elif [ ! -z "$left_reads" ] && [ ! -z "$right_reads" ] && [ "$quality_64" != 0 ] && [ "$tra_as" == 0 ] && [ "$tra_cuff" != 0 ]; then
    mkdir "$bam_out"
    numb=$(ls "${left_reads[@]}" | wc -l)
    for f in "${left_reads[@]}"; do
      extension=$(echo "$f" | sed -r 's/.*(fq|fq.gz|fastq|fastq.gz)$/\1/')
       
      if [[ "$extension" =~ "fq.gz" ]]; then
        filename=$(basename "$f" ".fq.gz")
        filename2=${filename/_R1/_R2}
        filename3=$(echo $filename | sed 's/_R1//')
        
        if [ $numb -eq 1 ] && [ "$cuff_merge" == 0 ]; then
            hisat2 -x $fbname --rna-strandness $lib_type -1 ${filename}.fq.gz -2 ${filename2}.fq.gz -S $filename3.sam -p $num_threads -5 $five_trim -3 $three_trim --min-intronlen $min_intl --max-intronlen $max_intl --dta --phred64
            cufflinks_non_SRA
            rm $fbname*
            exit
        elif [ $numb -eq 1 ] && [ "$cuff_merge" != 0 ]; then
            hisat2 -x $fbname --rna-strandness $lib_type -1 ${filename}.fq.gz -2 ${filename2}.fq.gz -S $filename3.sam -p $num_threads -5 $five_trim -3 $three_trim --min-intronlen $min_intl --max-intronlen $max_intl --dta --phred64
            cufflinks_non_SRA
            rm $fbname*
            echo "cuffmerge only works with more than one file"
            exit    
        elif [ $numb -gt 1 ] && [ "$cuff_merge" == 0 ]; then
            hisat2 -x $fbname --rna-strandness $lib_type -1 ${filename}.fq.gz -2 ${filename2}.fq.gz -S $filename3.sam -p $num_threads -5 $five_trim -3 $three_trim --min-intronlen $min_intl --max-intronlen $max_intl --dta --phred64
            cufflinks_non_SRA
        elif [ $numb -gt 1 ] && [ "$cuff_merge" != 0 ]; then
            hisat2 -x $fbname --rna-strandness $lib_type -1 ${filename}.fq.gz -2 ${filename2}.fq.gz -S $filename3.sam -p $num_threads -5 $five_trim -3 $three_trim --min-intronlen $min_intl --max-intronlen $max_intl --dta --phred64
            cufflinks_non_SRA
        fi

      elif [[ "$extension" =~ "fastq.gz" ]]; then
        filename=$(basename "$f" ".fastq.gz")
        filename2=${filename/_R1/_R2}
	      filename3=$(echo $filename | sed 's/_R1//')

        if [ $numb -eq 1 ] && [ "$cuff_merge" == 0 ]; then
            hisat2 -x $fbname --rna-strandness $lib_type -1 ${filename}.fastq.gz -2 ${filename2}.fastq.gz -S $filename3.sam -p $num_threads -5 $five_trim -3 $three_trim --min-intronlen $min_intl --max-intronlen $max_intl --dta --phred64
            cufflinks_non_SRA
            rm $fbname*
            exit
        elif [ $numb -eq 1 ] && [ "$cuff_merge" != 0 ]; then
            hisat2 -x $fbname --rna-strandness $lib_type -1 ${filename}.fastq.gz -2 ${filename2}.fastq.gz -S $filename3.sam -p $num_threads -5 $five_trim -3 $three_trim --min-intronlen $min_intl --max-intronlen $max_intl --dta --phred64
            cufflinks_non_SRA
            rm $fbname*
            echo "cuffmerge only works with more than one file"
            exit    
        elif [ $numb -gt 1 ] && [ "$cuff_merge" == 0 ]; then
            hisat2 -x $fbname --rna-strandness $lib_type -1 ${filename}.fastq.gz -2 ${filename2}.fastq.gz -S $filename3.sam -p $num_threads -5 $five_trim -3 $three_trim --min-intronlen $min_intl --max-intronlen $max_intl --dta --phred64
            cufflinks_non_SRA
        elif [ $numb -gt 1 ] && [ "$cuff_merge" != 0 ]; then
            hisat2 -x $fbname --rna-strandness $lib_type -1 ${filename}.fastq.gz -2 ${filename2}.fastq.gz -S $filename3.sam -p $num_threads -5 $five_trim -3 $three_trim --min-intronlen $min_intl --max-intronlen $max_intl --dta --phred64
            cufflinks_non_SRA
        fi

      elif [[ "$extension" =~ "fq" ]]; then
        filename=$(basename "$f" ".fq")
        filename2=${filename/_R1/_R2}
	      filename3=$(echo $filename | sed 's/_R1//')

        if [ $numb -eq 1 ] && [ "$cuff_merge" == 0 ]; then
            hisat2 -x $fbname --rna-strandness $lib_type -1 ${filename}.fq -2 ${filename2}.fq -S $filename3.sam -p $num_threads -5 $five_trim -3 $three_trim --min-intronlen $min_intl --max-intronlen $max_intl --dta --phred64
            cufflinks_non_SRA
            rm $fbname*
            exit
        elif [ $numb -eq 1 ] && [ "$cuff_merge" != 0 ]; then
            hisat2 -x $fbname --rna-strandness $lib_type -1 ${filename}.fq -2 ${filename2}.fq -S $filename3.sam -p $num_threads -5 $five_trim -3 $three_trim --min-intronlen $min_intl --max-intronlen $max_intl --dta --phred64
            cufflinks_non_SRA
            rm $fbname*
            echo "cuffmerge only works with more than one file"
            exit    
        elif [ $numb -gt 1 ] && [ "$cuff_merge" == 0 ]; then
            hisat2 -x $fbname --rna-strandness $lib_type -1 ${filename}.fq -2 ${filename2}.fq -S $filename3.sam -p $num_threads -5 $five_trim -3 $three_trim --min-intronlen $min_intl --max-intronlen $max_intl --dta --phred64
            cufflinks_non_SRA
        elif [ $numb -gt 1 ] && [ "$cuff_merge" != 0 ]; then
            hisat2 -x $fbname --rna-strandness $lib_type -1 ${filename}.fq -2 ${filename2}.fq -S $filename3.sam -p $num_threads -5 $five_trim -3 $three_trim --min-intronlen $min_intl --max-intronlen $max_intl --dta --phred64
            cufflinks_non_SRA
        fi

      elif [[ "$extension" =~ "fastq" ]]; then
        filename=$(basename "$f" ".fastq")
        filename2=${filename/_R1/_R2}
	      filename3=$(echo $filename | sed 's/_R1//')

        if [ $numb -eq 1 ] && [ "$cuff_merge" == 0 ]; then
            hisat2 -x $fbname --rna-strandness $lib_type -1 ${filename}.fastq -2 ${filename2}.fastq -S $filename3.sam -p $num_threads -5 $five_trim -3 $three_trim --min-intronlen $min_intl --max-intronlen $max_intl --dta --phred64
            cufflinks_non_SRA
            rm $fbname*
            exit
        elif [ $numb -eq 1 ] && [ "$cuff_merge" != 0 ]; then
            hisat2 -x $fbname --rna-strandness $lib_type -1 ${filename}.fastq -2 ${filename2}.fastq -S $filename3.sam -p $num_threads -5 $five_trim -3 $three_trim --min-intronlen $min_intl --max-intronlen $max_intl --dta --phred64
            cufflinks_non_SRA
            rm $fbname*
            echo "cuffmerge only works with more than one file"
            exit    
        elif [ $numb -gt 1 ] && [ "$cuff_merge" == 0 ]; then
            hisat2 -x $fbname --rna-strandness $lib_type -1 ${filename}.fastq -2 ${filename2}.fastq -S $filename3.sam -p $num_threads -5 $five_trim -3 $three_trim --min-intronlen $min_intl --max-intronlen $max_intl --dta --phred64
            cufflinks_non_SRA
        elif [ $numb -gt 1 ] && [ "$cuff_merge" != 0 ]; then
            hisat2 -x $fbname --rna-strandness $lib_type -1 ${filename}.fastq -2 ${filename2}.fastq -S $filename3.sam -p $num_threads -5 $five_trim -3 $three_trim --min-intronlen $min_intl --max-intronlen $max_intl --dta --phred64
            cufflinks_non_SRA
        fi   

      elif [ "$extension" != "fastq" ] || [ "$extension" != "fq" ] || [ "$extension" != "fastq.gz" ] || [ "$extension" != "fq.gz" ]; then
        rm $fbname*
        echo "The extension" "$extension" "is not supported. Only .fq, .fq.gz, .fastq, .fastq.gz are only supported" 1>&2        
        exit 64
      fi 
    done

    if [ $numb -gt 1 ] && [ "$cuff_merge" == 0 ]; then
      exit
    elif [ $numb -gt 1 ] && [ "$cuff_merge" != 0 ]; then
      cuff_merge_fun
    fi    


# ############################################################################################################################################################################################################################
# # Single end reads
# ############################################################################################################################################################################################################################

#### Stringtie

# Phred 33

elif [ ! -z "$single_reads" ] && [ "$quality_33" != 0 ] && [ "$tra_as" != 0 ] && [ "$tra_cuff" == 0 ]; then
    mkdir "$bam_out"
    numb=$(ls "${single_reads[@]}" | wc -l)
    for f in "${single_reads[@]}"; do
      if [ $numb -eq 1 ] && [ "$cuff_merge" == 0 ]; then
        extension=$(echo "$f" | sed -r 's/.*(fq|fq.gz|fastq|fastq.gz)$/\1/')
        filename=$(basename "$f" ".$extension")
        hisat2 -x $fbname --rna-strandness $lib_type -U $f -S $filename.sam -p $num_threads -5 $five_trim -3 $three_trim --min-intronlen $min_intl --max-intronlen $max_intl --dta --phred33
        stringtie_non_SRA_single
        rm $fbname*
        exit
      elif [ $numb -eq 1 ] && [ "$cuff_merge" != 0 ]; then
        extension=$(echo "$f" | sed -r 's/.*(fq|fq.gz|fastq|fastq.gz)$/\1/')
        filename=$(basename "$f" ".$extension")
        hisat2 -x $fbname --rna-strandness $lib_type -U $f -S $filename.sam -p $num_threads -5 $five_trim -3 $three_trim --min-intronlen $min_intl --max-intronlen $max_intl --dta --phred33
        stringtie_non_SRA_single
        rm $fbname*
        echo "cuffmerge only works with more than one file"
        exit
      elif [ $numb -gt 1 ] && [ "$cuff_merge" == 0 ]; then
        extension=$(echo "$f" | sed -r 's/.*(fq|fq.gz|fastq|fastq.gz)$/\1/')
        filename=$(basename "$f" ".$extension")
        hisat2 -x $fbname --rna-strandness $lib_type -U $f -S $filename.sam -p $num_threads -5 $five_trim -3 $three_trim --min-intronlen $min_intl --max-intronlen $max_intl --dta --phred33
        stringtie_non_SRA_single   
      elif [ $numb -gt 1 ] && [ "$cuff_merge" != 0 ]; then
        extension=$(echo "$f" | sed -r 's/.*(fq|fq.gz|fastq|fastq.gz)$/\1/')
        filename=$(basename "$f" ".$extension")
        hisat2 -x $fbname --rna-strandness $lib_type -U $f -S $filename.sam -p $num_threads -5 $five_trim -3 $three_trim --min-intronlen $min_intl --max-intronlen $max_intl --dta --phred33
        stringtie_non_SRA_single
      fi 
    done

    if [ $numb -gt 1 ] && [ "$cuff_merge" == 0 ]; then
      exit
    elif [ $numb -gt 1 ] && [ "$cuff_merge" != 0 ]; then
      cuff_merge_fun
    fi    

# Phred 64

elif [ ! -z "$single_reads" ] && [ "$quality_64" != 0 ] && [ "$tra_as" != 0 ] && [ "$tra_cuff" == 0 ]; then
    mkdir "$bam_out"
    numb=$(ls "${single_reads[@]}" | wc -l)
    for f in "${single_reads[@]}"; do
      if [ $numb -eq 1 ] && [ "$cuff_merge" == 0 ]; then
        extension=$(echo "$f" | sed -r 's/.*(fq|fq.gz|fastq|fastq.gz)$/\1/')
        filename=$(basename "$f" ".$extension")
        hisat2 -x $fbname --rna-strandness $lib_type -U $f -S $filename.sam -p $num_threads -5 $five_trim -3 $three_trim --min-intronlen $min_intl --max-intronlen $max_intl --dta --phred64
        stringtie_non_SRA_single
        rm $fbname*
        exit 
      elif [ $numb -eq 1 ] && [ "$cuff_merge" != 0 ]; then
        extension=$(echo "$f" | sed -r 's/.*(fq|fq.gz|fastq|fastq.gz)$/\1/')
        filename=$(basename "$f" ".$extension")
        hisat2 -x $fbname --rna-strandness $lib_type -U $f -S $filename.sam -p $num_threads -5 $five_trim -3 $three_trim --min-intronlen $min_intl --max-intronlen $max_intl --dta --phred64
        stringtie_non_SRA_single
        rm $fbname*
        echo "cuffmerge only works with more than one file"
        exit  
      elif [ $numb -gt 1 ] && [ "$cuff_merge" == 0 ]; then
        extension=$(echo "$f" | sed -r 's/.*(fq|fq.gz|fastq|fastq.gz)$/\1/')
        filename=$(basename "$f" ".$extension")
        hisat2 -x $fbname --rna-strandness $lib_type -U $f -S $filename.sam -p $num_threads -5 $five_trim -3 $three_trim --min-intronlen $min_intl --max-intronlen $max_intl --dta --phred64
        stringtie_non_SRA_single
      elif [ $numb -gt 1 ] && [ "$cuff_merge" != 0 ]; then
        extension=$(echo "$f" | sed -r 's/.*(fq|fq.gz|fastq|fastq.gz)$/\1/')
        filename=$(basename "$f" ".$extension")
        hisat2 -x $fbname --rna-strandness $lib_type -U $f -S $filename.sam -p $num_threads -5 $five_trim -3 $three_trim --min-intronlen $min_intl --max-intronlen $max_intl --dta --phred64
        stringtie_non_SRA_single
      fi 
    done
    
    if [ $numb -gt 1 ] && [ "$cuff_merge" == 0 ]; then
      exit
    elif [ $numb -gt 1 ] && [ "$cuff_merge" != 0 ]; then
      cuff_merge_fun
    fi 

# Cufflinks

# Phred 33

elif [ ! -z "$single_reads" ] && [ "$quality_33" != 0 ] && [ "$tra_as" == 0 ] && [ "$tra_cuff" != 0 ]; then
    mkdir "$bam_out"
    numb=$(ls "${single_reads[@]}" | wc -l)
    for f in "${single_reads[@]}"; do
      if [ $numb -eq 1 ] && [ "$cuff_merge" == 0 ]; then 
        extension=$(echo "$f" | sed -r 's/.*(fq|fq.gz|fastq|fastq.gz)$/\1/')
        filename=$(basename "$f" ".$extension")
        hisat2 -x $fbname --rna-strandness $lib_type -U $f -S $filename.sam -p $num_threads -5 $five_trim -3 $three_trim --min-intronlen $min_intl --max-intronlen $max_intl --dta-cufflinks --phred33
        cufflinks_non_SRA_single
        rm $fbname*
        exit
      elif [ $numb -eq 1 ] && [ "$cuff_merge" != 0 ]; then 
        extension=$(echo "$f" | sed -r 's/.*(fq|fq.gz|fastq|fastq.gz)$/\1/')
        filename=$(basename "$f" ".$extension")
        hisat2 -x $fbname --rna-strandness $lib_type -U $f -S $filename.sam -p $num_threads -5 $five_trim -3 $three_trim --min-intronlen $min_intl --max-intronlen $max_intl --dta-cufflinks --phred33
        cufflinks_non_SRA_single
        rm $fbname*
        echo "cuffmerge only works with more than one file"
        exit  
      elif [ $numb -gt 1 ] && [ "$cuff_merge" == 0 ]; then
        extension=$(echo "$f" | sed -r 's/.*(fq|fq.gz|fastq|fastq.gz)$/\1/')
        filename=$(basename "$f" ".$extension")
        hisat2 -x $fbname --rna-strandness $lib_type -U $f -S $filename.sam -p $num_threads -5 $five_trim -3 $three_trim --min-intronlen $min_intl --max-intronlen $max_intl --dta-cufflinks --phred33
        cufflinks_non_SRA_single
      elif [ $numb -gt 1 ] && [ "$cuff_merge" != 0 ]; then
        extension=$(echo "$f" | sed -r 's/.*(fq|fq.gz|fastq|fastq.gz)$/\1/')
        filename=$(basename "$f" ".$extension")
        hisat2 -x $fbname --rna-strandness $lib_type -U $f -S $filename.sam -p $num_threads -5 $five_trim -3 $three_trim --min-intronlen $min_intl --max-intronlen $max_intl --dta-cufflinks --phred33
        cufflinks_non_SRA_single
      fi          
    done
    
    if [ $numb -gt 1 ] && [ "$cuff_merge" == 0 ]; then
      exit
    elif [ $numb -gt 1 ] && [ "$cuff_merge" != 0 ]; then
      cuff_merge_fun
    fi

# Phred 33

elif [ ! -z "$single_reads" ] && [ "$quality_64" != 0 ] && [ "$tra_as" == 0 ] && [ "$tra_cuff" != 0 ]; then
    mkdir "$bam_out"
    numb=$(ls "${single_reads[@]}" | wc -l)
    for f in "${single_reads[@]}"; do
      if [ $numb -eq 1 ] && [ "$cuff_merge" == 0 ]; then 
        extension=$(echo "$f" | sed -r 's/.*(fq|fq.gz|fastq|fastq.gz)$/\1/')
        filename=$(basename "$f" ".$extension")
        hisat2 -x $fbname --rna-strandness $lib_type -U $f -S $filename.sam -p $num_threads -5 $five_trim -3 $three_trim --min-intronlen $min_intl --max-intronlen $max_intl --dta-cufflinks --phred64
        cufflinks_non_SRA_single
        rm $fbname*
        exit
      elif [ $numb -eq 1 ] && [ "$cuff_merge" != 0 ]; then 
        extension=$(echo "$f" | sed -r 's/.*(fq|fq.gz|fastq|fastq.gz)$/\1/')
        filename=$(basename "$f" ".$extension")
        hisat2 -x $fbname --rna-strandness $lib_type -U $f -S $filename.sam -p $num_threads -5 $five_trim -3 $three_trim --min-intronlen $min_intl --max-intronlen $max_intl --dta-cufflinks --phred64
        cufflinks_non_SRA_single
        rm $fbname*
        echo "cuffmerge only works with more than one file"
        exit  
      elif [ $numb -gt 1 ] && [ "$cuff_merge" == 0 ]; then
        extension=$(echo "$f" | sed -r 's/.*(fq|fq.gz|fastq|fastq.gz)$/\1/')
        filename=$(basename "$f" ".$extension")
        hisat2 -x $fbname --rna-strandness $lib_type -U $f -S $filename.sam -p $num_threads -5 $five_trim -3 $three_trim --min-intronlen $min_intl --max-intronlen $max_intl --dta-cufflinks --phred64
        cufflinks_non_SRA_single
      elif [ $numb -gt 1 ] && [ "$cuff_merge" != 0 ]; then
        extension=$(echo "$f" | sed -r 's/.*(fq|fq.gz|fastq|fastq.gz)$/\1/')
        filename=$(basename "$f" ".$extension")
        hisat2 -x $fbname --rna-strandness $lib_type -U $f -S $filename.sam -p $num_threads -5 $five_trim -3 $three_trim --min-intronlen $min_intl --max-intronlen $max_intl --dta-cufflinks --phred64
        cufflinks_non_SRA_single
      fi          
    done

    if [ $numb -gt 1 ] && [ "$cuff_merge" == 0 ]; then
      exit
    elif [ $numb -gt 1 ] && [ "$cuff_merge" != 0 ]; then
      cuff_merge_fun
    fi

# ############################################################################################################################################################################################################################
# # SRA reads
# ###########################################################################################################################################################################################################################

# Stringtie

# phred 33

elif [ ! -z $sra_id ] && [ "$quality_33" != 0 ] && [ "$tra_as" != 0 ] && [ "$tra_cuff" == 0 ]; then
  if [[ -f $sra_id ]]; then
    mkdir "$bam_out"    
    while read f; do
      hisat2 -x $fbname --rna-strandness $lib_type --sra-acc $f -S $f.sam -p $num_threads -5 $five_trim -3 $three_trim --min-intronlen $min_intl --max-intronlen $max_intl --dta --phred33
      stringtie_SRA_multi      
      done < "$sra_id" 
      if [ "$cuff_merge" == 0 ]; then
        exit
      elif [ "$cuff_merge" != 0 ]; then
        cuff_merge_fun
      fi
  else    
    mkdir "$bam_out"
    hisat2 -x $fbname --rna-strandness $lib_type --sra-acc $sra_id -S $sra_id.sam -p $num_threads -5 $five_trim -3 $three_trim --min-intronlen $min_intl --max-intronlen $max_intl --dta --phred33
    stringtie_SRA_single
    if [ "$cuff_merge" != 0 ]; then
      echo "cuffmerge only works with more than one SRA accesions. Use File containing SRA id's option to use cuffmerge"
    fi  
  fi

# phred 64

elif [ ! -z $sra_id ] && [ "$quality_64" != 0 ] && [ "$tra_as" != 0 ] && [ "$tra_cuff" == 0 ]; then
  if [[ -f $sra_id ]]; then
    mkdir "$bam_out"    
    while read f; do
      hisat2 -x $fbname --rna-strandness $lib_type --sra-acc $f -S $f.sam -p $num_threads -5 $five_trim -3 $three_trim --min-intronlen $min_intl --max-intronlen $max_intl --dta --phred64
      stringtie_SRA_multi      
      done < "$sra_id"
      if [ "$cuff_merge" == 0 ]; then
        exit
      elif [ "$cuff_merge" != 0 ]; then
        cuff_merge_fun
      fi
  else    
    mkdir "$bam_out"
    hisat2 -x $fbname --rna-strandness $lib_type --sra-acc $sra_id -S $sra_id.sam -p $num_threads -5 $five_trim -3 $three_trim --min-intronlen $min_intl --max-intronlen $max_intl --dta --phred64
    stringtie_SRA_single
    if [ "$cuff_merge" != 0 ]; then
      echo "cuffmerge only works with more than one SRA accesions. Use File containing SRA id's option to use cuffmerge"
    fi
  fi 

# Cufflinks

# phred 33

elif [ ! -z $sra_id ] && [ "$quality_33" != 0 ] && [ "$tra_as" == 0 ] && [ "$tra_cuff" != 0 ]; then
  if [[ -f $sra_id ]]; then
    mkdir "$bam_out"    
    while read f; do
      hisat2 -x $fbname --rna-strandness $lib_type --sra-acc $f -S $f.sam -p $num_threads -5 $five_trim -3 $three_trim --min-intronlen $min_intl --max-intronlen $max_intl --dta-cufflinks --phred33
      cufflinks_SRA_multi      
      done < "$sra_id"
      if [ "$cuff_merge" == 0 ]; then
        exit
      elif [ "$cuff_merge" != 0 ]; then
        cuff_merge_fun
      fi
  else    
    mkdir "$bam_out"
    hisat2 -x $fbname --rna-strandness $lib_type --sra-acc $sra_id -S $sra_id.sam -p $num_threads -5 $five_trim -3 $three_trim --min-intronlen $min_intl --max-intronlen $max_intl --dta-cufflinks --phred33
    cufflinks_SRA_single
    if [ "$cuff_merge" != 0 ]; then
      echo "cuffmerge only works with more than one SRA accesions. Use File containing SRA id's option to use cuffmerge"
    fi
  fi

# Phred 64

elif [ ! -z $sra_id ] && [ "$quality_64" != 0 ] && [ "$tra_as" == 0 ] && [ "$tra_cuff" != 0 ]; then
  if [[ -f $sra_id ]]; then
    mkdir "$bam_out"    
    while read f; do
      hisat2 -x $fbname --rna-strandness $lib_type --sra-acc $f -S $f.sam -p $num_threads -5 $five_trim -3 $three_trim --min-intronlen $min_intl --max-intronlen $max_intl --dta-cufflinks --phred64
      cufflinks_SRA_multi      
      done < "$sra_id"
      if [ "$cuff_merge" == 0 ]; then
        exit
      elif [ "$cuff_merge" != 0 ]; then
        cuff_merge_fun
      fi
  else    
    mkdir "$bam_out"
    hisat2 -x $fbname --rna-strandness $lib_type --sra-acc $sra_id -S $sra_id.sam -p $num_threads -5 $five_trim -3 $three_trim --min-intronlen $min_intl --max-intronlen $max_intl --dta-cufflinks --phred64
    cufflinks_SRA_single
    if [ "$cuff_merge" != 0 ]; then
      echo "cuffmerge only works with more than one SRA accesions. Use File containing SRA id's option to use cuffmerge"
    fi
  fi
fi

# Clean up the reference genomes
rm $fbname*
