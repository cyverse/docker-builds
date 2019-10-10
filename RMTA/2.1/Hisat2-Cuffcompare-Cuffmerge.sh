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


while getopts ":hg:i:A:l:1:2:U:O:s:p:5:3:f:qQtcem:M:k:" opt; do
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
    k)
    fthreshold=$OPTARG # FPKM filter that you would like to apply to identified transcripts
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
### No SRA #######
##################

### Paired end #####

# Coverage cut-off function for paired end reads

coverge_cuffoff_non_SRA()
{
    if [ "$threshold" -eq "$param5" ]; then
       grep " transcript" "$bam_out"/$filename3.gtf | grep -e 'cov "4\.' -e 'cov "3\.' -e 'cov "2\.' -e 'cov "1\.' -e 'cov "0\.' | cut -f 9 | cut -d " " -f 4 | sort -u > listtoremove.txt
       grep -vFf listtoremove.txt "$bam_out"/$filename3.gtf >"$bam_out"/$filename3.gtf.filtered.gtf
       var=$(grep -vFf listtoremove.txt "$bam_out"/$filename3.gtf.filtered.gtf | wc -l)
       if [ "$var" -eq 0 ]; then
          echo No transcripts have FPKM values exceeding your "$param5" cut-off have been found in $filename3.gtf.filtered.gtf. Try lowering your cut-off 1>&2
       fi
    elif [ "$threshold" -eq "$param4" ]; then
       grep " transcript" "$bam_out"/$filename3.gtf | grep -e 'cov "3\.' -e 'cov "2\.' -e 'cov "1\.' -e 'cov "0\.' | cut -f 9 | cut -d " " -f 4 | sort -u > listtoremove.txt
       grep -vFf listtoremove.txt "$bam_out"/$filename3.gtf >"$bam_out"/$filename3.gtf.filtered.gtf
       var=$(grep -vFf listtoremove.txt "$bam_out"/$filename3.gtf.filtered.gtf | wc -l)
       if [ "$var" -eq 0 ]; then
          echo No transcripts have FPKM values exceeding your "$param4" cut-off have been found in $filename3.gtf.filtered.gtf. Try lowering your cut-off 1>&2
       fi
    elif [ "$threshold" -eq "$param3" ]; then
       grep " transcript" "$bam_out"/$filename3.gtf | grep -e 'cov "2\.' -e 'cov "1\.' -e 'cov "0\.' | cut -f 9 | cut -d " " -f 4 | sort -u > listtoremove.txt
       grep -vFf listtoremove.txt "$bam_out"/$filename3.gtf >"$bam_out"/$filename3.gtf.filtered.gtf
       var=$(grep -vFf listtoremove.txt "$bam_out"/$filename3.gtf.filtered.gtf | wc -l)
       if [ "$var" -eq 0 ]; then
          echo No transcripts have FPKM values exceeding your "$param3" cut-off have been found in $filename3.gtf.filtered.gtf. Try lowering your cut-off 1>&2
       fi
    elif [ "$threshold" -eq "$param2" ]; then
       grep " transcript" "$bam_out"/$filename3.gtf | grep -e 'cov "1\.' -e 'cov "0\.' | cut -f 9 | cut -d " " -f 4 | sort -u > listtoremove.txt
       grep -vFf listtoremove.txt "$bam_out"/$filename3.gtf >"$bam_out"/$filename3.gtf.filtered.gtf
       var=$(grep -vFf listtoremove.txt "$bam_out"/$filename3.gtf.filtered.gtf | wc -l)
       if [ "$var" -eq 0 ]; then
          echo No transcripts have FPKM values exceeding your "$param2" cut-off have been found in $filename3.gtf.filtered.gtf. Try lowering your cut-off 1>&2
       fi
    elif [ "$threshold" -eq "$param1" ]; then
       grep " transcript" "$bam_out"/$filename3.gtf | grep -e 'cov "0\.' | cut -f 9 | cut -d " " -f 4 | sort -u > listtoremove.txt
       grep -vFf listtoremove.txt "$bam_out"/$filename3.gtf >"$bam_out"/$filename3.gtf.filtered.gtf
       var=$(grep -vFf listtoremove.txt "$bam_out"/$filename3.gtf.filtered.gtf | wc -l)
       if [ "$var" -eq 0 ]; then
          echo No transcripts have FPKM values exceeding your "$param1" cut-off have been found in $filename3.gtf.filtered.gtf. Try lowering your cut-off 1>&2
       fi
    elif [ "$threshold" -eq "$param0" ]; then
       grep " transcript" "$bam_out"/$filename3.gtf | grep -e 'cov "0\.000' | cut -f 9 | cut -d " " -f 4 | sort -u > listtoremove.txt
       grep -vFf listtoremove.txt "$bam_out"/$filename3.gtf >"$bam_out"/$filename3.gtf.filtered.gtf
    else
        echo "Invalid coverage parameter. Please select a whole number between 0-5"
        exit
    fi
}

# FPKM cut-off function for paired end reads

fpkm_cuffoff_non_SRA()
{
    if [ "$fthreshold" -eq "$param5" ]; then
       grep " transcript" "$bam_out"/$filename3.gtf | grep -e 'FPKM "4\.' -e 'FPKM "3\.' -e 'FPKM "2\.' -e 'FPKM "1\.' -e 'FPKM "0\.' | cut -f 9 | cut -d " " -f 4 | sort -u > listtoremove.txt
       grep -vFf listtoremove.txt "$bam_out"/$filename3.gtf >"$bam_out"/$filename3.gtf.filtered.gtf
       var=$(grep -vFf listtoremove.txt "$bam_out"/$filename3.gtf.filtered.gtf | wc -l)
       if [ "$var" -eq 0 ]; then
          echo No transcripts have FPKM values exceeding your "$param5" cut-off have been found in $filename3.gtf.filtered.gtf. Try lowering your cut-off 1>&2
       fi
    elif [ "$fthreshold" -eq "$param4" ]; then
       grep " transcript" "$bam_out"/$filename3.gtf | grep -e 'FPKM "3\.' -e 'FPKM "2\.' -e 'FPKM "1\.' -e 'FPKM "0\.' | cut -f 9 | cut -d " " -f 4 | sort -u > listtoremove.txt
       grep -vFf listtoremove.txt "$bam_out"/$filename3.gtf >"$bam_out"/$filename3.gtf.filtered.gtf
       var=$(grep -vFf listtoremove.txt "$bam_out"/$filename3.gtf.filtered.gtf | wc -l)
       if [ "$var" -eq 0 ]; then
          echo No transcripts have FPKM values exceeding your "$param4" cut-off have been found in $filename3.gtf.filtered.gtf. Try lowering your cut-off 1>&2
       fi
    elif [ "$fthreshold" -eq "$param3" ]; then
       grep " transcript" "$bam_out"/$filename3.gtf | grep -e 'FPKM "2\.' -e 'FPKM "1\.' -e 'FPKM "0\.' | cut -f 9 | cut -d " " -f 4 | sort -u > listtoremove.txt
       grep -vFf listtoremove.txt "$bam_out"/$filename3.gtf >"$bam_out"/$filename3.gtf.filtered.gtf
       var=$(grep -vFf listtoremove.txt "$bam_out"/$filename3.gtf.filtered.gtf | wc -l)
       if [ "$var" -eq 0 ]; then
          echo No transcripts have FPKM values exceeding your "$param3" cut-off have been found in $filename3.gtf.filtered.gtf. Try lowering your cut-off 1>&2
       fi
    elif [ "$fthreshold" -eq "$param2" ]; then
       grep " transcript" "$bam_out"/$filename3.gtf | grep -e 'FPKM "1\.' -e 'FPKM "0\.' | cut -f 9 | cut -d " " -f 4 | sort -u > listtoremove.txt
       grep -vFf listtoremove.txt "$bam_out"/$filename3.gtf >"$bam_out"/$filename3.gtf.filtered.gtf
       var=$(grep -vFf listtoremove.txt "$bam_out"/$filename3.gtf.filtered.gtf | wc -l)
       if [ "$var" -eq 0 ]; then
          echo No transcripts have FPKM values exceeding your "$param2" cut-off have been found in $filename3.gtf.filtered.gtf. Try lowering your cut-off 1>&2
       fi
    elif [ "$fthreshold" -eq "$param1" ]; then
       grep " transcript" "$bam_out"/$filename3.gtf | grep -e 'FPKM "0\.' | cut -f 9 | cut -d " " -f 4 | sort -u > listtoremove.txt
       grep -vFf listtoremove.txt "$bam_out"/$filename3.gtf >"$bam_out"/$filename3.gtf.filtered.gtf
       var=$(grep -vFf listtoremove.txt "$bam_out"/$filename3.gtf.filtered.gtf | wc -l)
       if [ "$var" -eq 0 ]; then
          echo No transcripts have FPKM values exceeding your "$param1" cut-off have been found in $filename3.gtf.filtered.gtf. Try lowering your cut-off 1>&2
       fi
    elif [ "$fthreshold" -eq "$param0" ]; then
       grep " transcript" "$bam_out"/$filename3.gtf | grep -e 'FPKM "0\.000' | cut -f 9 | cut -d " " -f 4 | sort -u > listtoremove.txt
       grep -vFf listtoremove.txt "$bam_out"/$filename3.gtf >"$bam_out"/$filename3.gtf.filtered.gtf
    else
        echo "Invalid coverage parameter. Please select a whole number between 0-5"
        exit
    fi
}

# Stringtie function for paired end reads

stringtie_non_SRA() 
{
    echo "###################"
    echo "Running Stringtie"
    echo "####################"
    echo "sambamba-0.6.8-linux-static sort --tmpdir=temp -t $num_threads -o "$bam_out"/$filename3.sorted.bam $filename3.bam"
    sambamba-0.6.8-linux-static sort --tmpdir=temp -t $num_threads -o "$bam_out"/$filename3.sorted.bam $filename3.bam
    rm -r temp $filename3.bam
    echo "stringtie -G $referenceannotation "$bam_out"/$filename3.sorted.bam -o "$bam_out"/$filename3.gtf -p $num_threads"
    stringtie -G $referenceannotation "$bam_out"/$filename3.sorted.bam -o "$bam_out"/$filename3.gtf -p $num_threads
    coverge_cuffoff_non_SRA
    fpkm_cuffoff_non_SRA 
    var=$(grep -vFf listtoremove.txt "$bam_out"/$filename3.gtf.filtered.gtf | wc -l)
       if [ "$var" -eq 0 ]; then
          rm listtoremove.txt "$bam_out"/$filename3.gtf.filtered.gtf
        else
          echo "###################"
          echo "Running Cuffcompare"
          echo "####################"
          echo "cuffcompare "$bam_out"/$filename3.gtf.filtered.gtf -r $referenceannotation -o $filename3"
          cuffcompare "$bam_out"/$filename3.gtf.filtered.gtf -r $referenceannotation -o $filename3
          rm listtoremove.txt
          mv *.tracking *.loci *.combined.gtf *.stats "$bam_out"
       fi
}

# Cufflinks function for paired end reads

cufflinks_non_SRA()
{
    echo "###################"
    echo "Running Cufflinks"
    echo "####################"
    echo "sambamba-0.6.8-linux-static sort --tmpdir=temp -t $num_threads -o "$bam_out"/$filename3.sorted.bam $filename3.bam"
    sambamba-0.6.8-linux-static sort --tmpdir=temp -t $num_threads -o "$bam_out"/$filename3.sorted.bam $filename3.bam
    rm -r temp $filename3.bam
    echo "cufflinks "$bam_out"/$filename3.sorted.bam -p $num_threads -g $referenceannotation -o "$bam_out""
    cufflinks "$bam_out"/$filename3.sorted.bam -p $num_threads -g $referenceannotation -o "$bam_out"
    mv "$bam_out"/skipped.gtf "$bam_out"/$filename3.skipped.gtf
    mv "$bam_out"/transcripts.gtf "$bam_out"/$filename3.gtf
    mv "$bam_out"/isoforms.fpkm_tracking "$bam_out"/$filename3.isoforms.fpkm_tracking
    mv "$bam_out"/genes.fpkm_tracking "$bam_out"/$filename3.genes.fpkm_tracking
    coverge_cuffoff_non_SRA
    fpkm_cuffoff_non_SRA
    var=$(grep -vFf listtoremove.txt "$bam_out"/$filename3.gtf.filtered.gtf | wc -l)
       if [ "$var" -eq 0 ]; then
          rm listtoremove.txt "$bam_out"/$filename3.gtf.filtered.gtf
        else
          echo "###################"
          echo "Running Cuffcompare"
          echo "####################"
          cuffcompare "$bam_out"/$filename3.gtf.filtered.gtf -r $referenceannotation -o $filename3
          rm listtoremove.txt
          mv *.tracking *.loci *.combined.gtf *.stats "$bam_out"
       fi  
}

### single end ####

# Coverage cut-off function for single end reads

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
 
# FPKM cut-off function for paired end reads

fpkm_cuffoff_non_SRA_single()
{
    if [ "$fthreshold" -eq "$param5" ]; then
       grep " transcript" "$bam_out"/$filename.gtf | grep -e 'FPKM "4.' -e 'FPKM "3.' -e 'FPKM "2.' -e 'FPKM "1.' -e 'FPKM "0.' | cut -f 9 | cut -d " " -f 4 | sort -u > listtoremove.txt
       grep -vFf listtoremove.txt "$bam_out"/$filename.gtf >"$bam_out"/$filename.gtf.filtered.gtf
       var=$(grep -vFf listtoremove.txt "$bam_out"/$filename.gtf.filtered.gtf | wc -l)
       if [ "$var" -eq 0 ]; then
          echo No transcripts have FPKM values exceeding your "$param5" cut-off have been found in $filename.gtf.filtered.gtf. Try lowering your cut-off 1>&2
       fi       
    elif [ "$fthreshold" -eq "$param4" ]; then
       grep " transcript" "$bam_out"/$filename.gtf | grep -e 'FPKM "3.' -e 'FPKM "2.' -e 'FPKM "1.' -e 'FPKM "0.' | cut -f 9 | cut -d " " -f 4 | sort -u > listtoremove.txt
       grep -vFf listtoremove.txt "$bam_out"/$filename.gtf >"$bam_out"/$filename.gtf.filtered.gtf
       var=$(grep -vFf listtoremove.txt "$bam_out"/$filename.gtf.filtered.gtf | wc -l)
       if [ "$var" -eq 0 ]; then
          echo No transcripts have FPKM values exceeding your "$param4" cut-off have been found in $filename.gtf.filtered.gtf. Try lowering your cut-off 1>&2
       fi
    elif [ "$fthreshold" -eq "$param3" ]; then
       grep " transcript" "$bam_out"/$filename.gtf | grep -e 'FPKM "2.' -e 'FPKM "1.' -e 'FPKM "0.' | cut -f 9 | cut -d " " -f 4 | sort -u > listtoremove.txt
       grep -vFf listtoremove.txt "$bam_out"/$filename.gtf >"$bam_out"/$filename.gtf.filtered.gtf
       var=$(grep -vFf listtoremove.txt "$bam_out"/$filename.gtf.filtered.gtf | wc -l)
       if [ "$var" -eq 0 ]; then
          echo No transcripts have FPKM values exceeding your "$param3" cut-off have been found in $filename.gtf.filtered.gtf. Try lowering your cut-off 1>&2
       fi
    elif [ "$fthreshold" -eq "$param2" ]; then
       grep " transcript" "$bam_out"/$filename.gtf | grep -e 'FPKM "1.' -e 'FPKM "0.' | cut -f 9 | cut -d " " -f 4 | sort -u > listtoremove.txt
       grep -vFf listtoremove.txt "$bam_out"/$filename.gtf >"$bam_out"/$filename.gtf.filtered.gtf
       var=$(grep -vFf listtoremove.txt "$bam_out"/$filename.gtf.filtered.gtf | wc -l)
       if [ "$var" -eq 0 ]; then
          echo No transcripts have FPKM values exceeding your "$param2" cut-off have been found in $filename.gtf.filtered.gtf. Try lowering your cut-off 1>&2
       fi
    elif [ "$fthreshold" -eq "$param1" ]; then
       grep " transcript" "$bam_out"/$filename.gtf | grep -e 'FPKM "0.' | cut -f 9 | cut -d " " -f 4 | sort -u > listtoremove.txt
       grep -vFf listtoremove.txt "$bam_out"/$filename.gtf >"$bam_out"/$filename.gtf.filtered.gtf
       var=$(grep -vFf listtoremove.txt "$bam_out"/$filename.gtf.filtered.gtf | wc -l)
       if [ "$var" -eq 0 ]; then
          echo No transcripts have FPKM values exceeding your "$param1" cut-off have been found in $filename.gtf.filtered.gtf. Try lowering your cut-off 1>&2
       fi
    elif [ "$fthreshold" -eq "$param0" ]; then
       grep " transcript" "$bam_out"/$filename.gtf | grep -e 'FPKM "0.000' | cut -f 9 | cut -d " " -f 4 | sort -u > listtoremove.txt
       grep -vFf listtoremove.txt "$bam_out"/$filename.gtf >"$bam_out"/$filename.gtf.filtered.gtf
    else
        echo "Invalid coverage parameter. Please select a whole number between 0-5"
        exit
    fi
}

# Stringtie function for paired end reads

stringtie_non_SRA_single() 
{
    echo "###################"
    echo "Running Stringtie"
    echo "####################"
    echo "sambamba-0.6.8-linux-static sort --tmpdir=temp -t $num_threads -o "$bam_out"/$filename.sorted.bam $filename.bam"
    sambamba-0.6.8-linux-static sort --tmpdir=temp -t $num_threads -o "$bam_out"/$filename.sorted.bam $filename.bam
    rm -r temp $filename.bam
    echo "stringtie -G $referenceannotation "$bam_out"/$filename.sorted.bam -o "$bam_out"/$filename.gtf -p $num_threads"
    stringtie -G $referenceannotation "$bam_out"/$filename.sorted.bam -o "$bam_out"/$filename.gtf -p $num_threads
    coverge_cuffoff_non_SRA_single
    fpkm_cuffoff_non_SRA_single
    var=$(grep -vFf listtoremove.txt "$bam_out"/$filename.gtf.filtered.gtf | wc -l)
       if [ "$var" -eq 0 ]; then
          rm listtoremove.txt "$bam_out"/$filename.gtf.filtered.gtf
        else
          echo "###################"
          echo "Running Cuffcompare"
          echo "####################"
          echo "cuffcompare "$bam_out"/$filename.gtf.filtered.gtf -r $referenceannotation -o $filename"
          cuffcompare "$bam_out"/$filename.gtf.filtered.gtf -r $referenceannotation -o $filename
          rm listtoremove.txt
          mv *.tracking *.loci *.combined.gtf *.stats "$bam_out"
       fi
}

# Cufflinks function for paired end reads

cufflinks_non_SRA_single()
{
    echo "###################"
    echo "Running Cufflinks"
    echo "####################"
    echo "sambamba-0.6.8-linux-static sort --tmpdir=temp -t $num_threads -o "$bam_out"/$filename.sorted.bam $filename.bam"
    sambamba-0.6.8-linux-static sort --tmpdir=temp -t $num_threads -o "$bam_out"/$filename.sorted.bam $filename.bam
    rm -r temp $filename.bam
    echo "cufflinks "$bam_out"/$filename.sorted.bam -p $num_threads -g $referenceannotation -o "$bam_out""
    cufflinks "$bam_out"/$filename.sorted.bam -p $num_threads -g $referenceannotation -o "$bam_out"
    mv "$bam_out"/skipped.gtf "$bam_out"/$filename.skipped.gtf
    mv "$bam_out"/transcripts.gtf "$bam_out"/$filename.gtf
    mv "$bam_out"/isoforms.fpkm_tracking "$bam_out"/$filename.isoforms.fpkm_tracking
    mv "$bam_out"/genes.fpkm_tracking "$bam_out"/$filename.genes.fpkm_tracking
    coverge_cuffoff_non_SRA_single
    fpkm_cuffoff_non_SRA_single        
    var=$(grep -vFf listtoremove.txt "$bam_out"/$filename.gtf.filtered.gtf | wc -l)
       if [ "$var" -eq 0 ]; then
          rm listtoremove.txt "$bam_out"/$filename.gtf.filtered.gtf
        else
          echo "###################"
          echo "Running Cuffcompare"
          echo "####################"
          echo "cuffcompare "$bam_out"/$filename.gtf.filtered.gtf -r $referenceannotation -o $filename"
          cuffcompare "$bam_out"/$filename.gtf.filtered.gtf -r $referenceannotation -o $filename
          rm listtoremove.txt
          mv *.tracking *.loci *.combined.gtf *.stats "$bam_out"
       fi
}

##################
### Single SRA ###
##################

# Coverage cut-off function for single SRAs

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

# FPKM cut-off cunction for Single SRAs

fpkm_cuffoff_SRA_single()
{
    if [ "$fthreshold" -eq "$param5" ]; then
       grep " transcript" "$bam_out"/$sra_id.gtf | grep -e 'FPKM "4.' -e 'FPKM "3.' -e 'FPKM "2.' -e 'FPKM "1.' -e 'FPKM "0.' | cut -f 9 | cut -d " " -f 4 | sort -u > listtoremove.txt
       grep -vFf listtoremove.txt "$bam_out"/$sra_id.gtf >"$bam_out"/$sra_id.gtf.filtered.gtf
       var=$(grep -vFf listtoremove.txt "$bam_out"/$sra_id.gtf.filtered.gtf | wc -l)
       if [ "$var" -eq 0 ]; then
          echo No transcripts have FPKM values exceeding your "$param5" cut-off have been found in $sra_id.gtf.filtered.gtf. Try lowering your cut-off 1>&2
       fi             
    elif [ "$fthreshold" -eq "$param4" ]; then
       grep " transcript" "$bam_out"/$sra_id.gtf | grep -e 'FPKM "3.' -e 'FPKM "2.' -e 'FPKM "1.' -e 'FPKM "0.' | cut -f 9 | cut -d " " -f 4 | sort -u > listtoremove.txt
       grep -vFf listtoremove.txt "$bam_out"/$sra_id.gtf >"$bam_out"/$sra_id.gtf.filtered.gtf
       var=$(grep -vFf listtoremove.txt "$bam_out"/$sra_id.gtf.filtered.gtf | wc -l)
       if [ "$var" -eq 0 ]; then
          echo No transcripts have FPKM values exceeding your "$param4" cut-off have been found in $sra_id.gtf.filtered.gtf. Try lowering your cut-off 1>&2
       fi             
    elif [ "$fthreshold" -eq "$param3" ]; then
       grep " transcript" "$bam_out"/$sra_id.gtf | grep -e 'FPKM "2.' -e 'FPKM "1.' -e 'FPKM "0.' | cut -f 9 | cut -d " " -f 4 | sort -u > listtoremove.txt
       grep -vFf listtoremove.txt "$bam_out"/$sra_id.gtf >"$bam_out"/$sra_id.gtf.filtered.gtf
       var=$(grep -vFf listtoremove.txt "$bam_out"/$sra_id.gtf.filtered.gtf | wc -l)
       if [ "$var" -eq 0 ]; then
          echo No transcripts have FPKM values exceeding your "$param3" cut-off have been found in $sra_id.gtf.filtered.gtf. Try lowering your cut-off 1>&2
       fi             
    elif [ "$fthreshold" -eq "$param2" ]; then
       grep " transcript" "$bam_out"/$sra_id.gtf | grep -e 'FPKM "1.' -e 'FPKM "0.' | cut -f 9 | cut -d " " -f 4 | sort -u > listtoremove.txt
       grep -vFf listtoremove.txt "$bam_out"/$sra_id.gtf >"$bam_out"/$sra_id.gtf.filtered.gtf
       var=$(grep -vFf listtoremove.txt "$bam_out"/$sra_id.gtf.filtered.gtf | wc -l)
       if [ "$var" -eq 0 ]; then
          echo No transcripts have FPKM values exceeding your "$param2" cut-off have been found in $sra_id.gtf.filtered.gtf. Try lowering your cut-off 1>&2
       fi             
    elif [ "$fthreshold" -eq "$param1" ]; then
       grep " transcript" "$bam_out"/$sra_id.gtf | grep -e 'FPKM "0.' | cut -f 9 | cut -d " " -f 4 | sort -u > listtoremove.txt
       grep -vFf listtoremove.txt "$bam_out"/$sra_id.gtf >"$bam_out"/$sra_id.gtf.filtered.gtf
       var=$(grep -vFf listtoremove.txt "$bam_out"/$sra_id.gtf.filtered.gtf | wc -l)
       if [ "$var" -eq 0 ]; then
          echo No transcripts have FPKM values exceeding your "$param1" cut-off have been found in $sra_id.gtf.filtered.gtf. Try lowering your cut-off 1>&2
       fi             
    elif [ "$fthreshold" -eq "$param0" ]; then
       grep " transcript" "$bam_out"/$sra_id.gtf | grep -e 'FPKM "0.000' | cut -f 9 | cut -d " " -f 4 | sort -u > listtoremove.txt
       grep -vFf listtoremove.txt "$bam_out"/$sra_id.gtf >"$bam_out"/$sra_id.gtf.filtered.gtf
    else
        echo "Invalid coverage parameter. Please select a whole number between 0-5"
        exit
    fi
}

# Stringtie function for single SRAs

stringtie_SRA_single()
{
    echo "###################"
    echo "Running Stringtie"
    echo "####################"
    echo "sambamba-0.6.8-linux-static sort --tmpdir=temp -t $num_threads -o "$bam_out"/$sra_id.sorted.bam $sra_id.bam"
    sambamba-0.6.8-linux-static sort --tmpdir=temp -t $num_threads -o "$bam_out"/$sra_id.sorted.bam $sra_id.bam
    rm -r temp $sra_id.bam
    echo "stringtie -G $referenceannotation "$bam_out"/$sra_id.sorted.bam -o "$bam_out"/$sra_id.gtf -p $num_threads"
    stringtie -G $referenceannotation "$bam_out"/$sra_id.sorted.bam -o "$bam_out"/$sra_id.gtf -p $num_threads
    coverge_cuffoff_SRA_single
    fpkm_cuffoff_SRA_single  
    var=$(grep -vFf listtoremove.txt "$bam_out"/$sra_id.gtf.filtered.gtf | wc -l)
    if [ "$var" -eq 0 ]; then
      rm listtoremove.txt "$bam_out"/$sra_id.gtf.filtered.gtf
    else
      echo "###################"
      echo "Running Cuffcompare"
      echo "####################"
      echo "cuffcompare "$bam_out"/$sra_id.gtf.filtered.gtf -r $referenceannotation -o $sra_id"
      cuffcompare "$bam_out"/$sra_id.gtf.filtered.gtf -r $referenceannotation -o $sra_id
      rm listtoremove.txt
      mv *.tracking *.loci *.combined.gtf *.stats "$bam_out"
   fi
}

# Cufflinks function for single SRAs

cufflinks_SRA_single()
{
    echo "###################"
    echo "Running Cufflinks"
    echo "####################"
    echo "sambamba-0.6.8-linux-static sort --tmpdir=temp -t $num_threads -o "$bam_out"/$sra_id.sorted.bam $sra_id.bam"
    sambamba-0.6.8-linux-static sort --tmpdir=temp -t $num_threads -o "$bam_out"/$sra_id.sorted.bam $sra_id.bam
    rm -r temp $sra_id.bam
    echo "cufflinks "$bam_out"/$sra_id.sorted.bam -p $num_threads -g $referenceannotation -o "$bam_out""
    cufflinks "$bam_out"/$sra_id.sorted.bam -p $num_threads -g $referenceannotation -o "$bam_out"
    mv "$bam_out"/transcripts.gtf "$bam_out"/$sra_id.gtf
    mv "$bam_out"/isoforms.fpkm_tracking "$bam_out"/$sra_id.isoforms.fpkm_tracking
    mv "$bam_out"/genes.fpkm_tracking "$bam_out"/$sra_id.genes.fpkm_tracking
    coverge_cuffoff_SRA_single
    fpkm_cuffoff_SRA_single           
    if [ "$var" -eq 0 ]; then
      rm listtoremove.txt "$bam_out"/$sra_id.gtf.filtered.gtf
    else
      echo "###################"
      echo "Running Cuffcompare"
      echo "####################"
      echo "cuffcompare "$bam_out"/$sra_id.gtf.filtered.gtf -r $referenceannotation -o $sra_id"
      cuffcompare "$bam_out"/$sra_id.gtf.filtered.gtf -r $referenceannotation -o $sra_id
      rm listtoremove.txt
      mv *.tracking *.loci *.combined.gtf *.stats "$bam_out"
   fi

}

##################
### Multi SRA ###
##################

# Coverage cut-off function for multiple SRAs

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

# FPKM cut-off function for multiple SRAs

fpkm_cuffoff_SRA_multi()
{
    if [ "$fthreshold" -eq "$param5" ]; then
       grep " transcript" "$bam_out"/$f.gtf | grep -e 'FPKM "4.' -e 'FPKM "3.' -e 'FPKM "2.' -e 'FPKM "1.' -e 'FPKM "0.' | cut -f 9 | cut -d " " -f 4 | sort -u > listtoremove.txt
       grep -vFf listtoremove.txt "$bam_out"/$f.gtf >"$bam_out"/$f.gtf.filtered.gtf
       var=$(grep -vFf listtoremove.txt "$bam_out"/$f.gtf.filtered.gtf | wc -l)
       if [ "$var" -eq 0 ]; then
          echo No transcripts have FPKM values exceeding your "$param5" cut-off have been found in $f.gtf.filtered.gtf. Try lowering your cut-off 1>&2
       fi      
    elif [ "$fthreshold" -eq "$param4" ]; then
       grep " transcript" "$bam_out"/$f.gtf | grep -e 'FPKM "3.' -e 'FPKM "2.' -e 'FPKM "1.' -e 'FPKM "0.' | cut -f 9 | cut -d " " -f 4 | sort -u > listtoremove.txt
       grep -vFf listtoremove.txt "$bam_out"/$f.gtf >"$bam_out"/$f.gtf.filtered.gtf
       var=$(grep -vFf listtoremove.txt "$bam_out"/$f.gtf.filtered.gtf | wc -l)
       if [ "$var" -eq 0 ]; then
          echo No transcripts have FPKM values exceeding your "$param4" cut-off have been found in $f.gtf.filtered.gtf. Try lowering your cut-off 1>&2
       fi       
    elif [ "$fthreshold" -eq "$param3" ]; then
       grep " transcript" "$bam_out"/$f.gtf | grep -e 'FPKM "2.' -e 'FPKM "1.' -e 'FPKM "0.' | cut -f 9 | cut -d " " -f 4 | sort -u > listtoremove.txt
       grep -vFf listtoremove.txt "$bam_out"/$f.gtf >"$bam_out"/$f.gtf.filtered.gtf
       var=$(grep -vFf listtoremove.txt "$bam_out"/$f.gtf.filtered.gtf | wc -l)
       if [ "$var" -eq 0 ]; then
          echo No transcripts have FPKM values exceeding your "$param3" cut-off have been found in $f.gtf.filtered.gtf. Try lowering your cut-off 1>&2
       fi      
    elif [ "$fthreshold" -eq "$param2" ]; then
       grep " transcript" "$bam_out"/$f.gtf | grep -e 'FPKM "1.' -e 'FPKM "0.' | cut -f 9 | cut -d " " -f 4 | sort -u > listtoremove.txt
       grep -vFf listtoremove.txt "$bam_out"/$f.gtf >"$bam_out"/$f.gtf.filtered.gtf
       var=$(grep -vFf listtoremove.txt "$bam_out"/$f.gtf.filtered.gtf | wc -l)
       if [ "$var" -eq 0 ]; then
          echo No transcripts have FPKM values exceeding your "$param2" cut-off have been found in $f.gtf.filtered.gtf. Try lowering your cut-off 1>&2
       fi       
    elif [ "$fthreshold" -eq "$param1" ]; then
       grep " transcript" "$bam_out"/$f.gtf | grep -e 'FPKM "0.' | cut -f 9 | cut -d " " -f 4 | sort -u > listtoremove.txt
       grep -vFf listtoremove.txt "$bam_out"/$f.gtf >"$bam_out"/$f.gtf.filtered.gtf
       var=$(grep -vFf listtoremove.txt "$bam_out"/$f.gtf.filtered.gtf | wc -l)
       if [ "$var" -eq 0 ]; then
          echo No transcripts have FPKM values exceeding your "$param1" cut-off have been found in $f.gtf.filtered.gtf. Try lowering your cut-off 1>&2
       fi       
    elif [ "$fthreshold" -eq "$param0" ]; then
       grep " transcript" "$bam_out"/$f.gtf | grep -e 'cov "0.000' | cut -f 9 | cut -d " " -f 4 | sort -u > listtoremove.txt
       grep -vFf listtoremove.txt "$bam_out"/$f.gtf >"$bam_out"/$f.gtf.filtered.gtf
    else
        echo "Invalid coverage parameter. Please select a whole number between 0-5"
        exit
    fi
}

# Stringtie function for multiple SRAs

stringtie_SRA_multi()
{
    echo "###################"
    echo "Running Stringtie"
    echo "####################"
    echo "sambamba-0.6.8-linux-static sort --tmpdir=temp -t $num_threads -o "$bam_out"/$f.sorted.bam $f.bam"
    sambamba-0.6.8-linux-static sort --tmpdir=temp -t $num_threads -o "$bam_out"/$f.sorted.bam $f.bam
    rm -r temp $f.bam
    echo "stringtie -G $referenceannotation "$bam_out"/$f.sorted.bam -o "$bam_out"/$f.gtf -p $num_threads"
    stringtie -G $referenceannotation "$bam_out"/$f.sorted.bam -o "$bam_out"/$f.gtf -p $num_threads
    coverge_cuffoff_SRA_multi
    fpkm_cuffoff_SRA_multi
    var=$(grep -vFf listtoremove.txt "$bam_out"/$f.gtf.filtered.gtf | wc -l)
    if [ "$var" -eq 0 ]; then
      rm listtoremove.txt "$bam_out"/$f.gtf.filtered.gtf
    else
      echo "###################"
      echo "Running Cuffcompare"
      echo "####################"
      echo "cuffcompare "$bam_out"/$f.gtf.filtered.gtf -r $referenceannotation -o $f"
      cuffcompare "$bam_out"/$f.gtf.filtered.gtf -r $referenceannotation -o $f
      rm listtoremove.txt
      mv *.tracking *.loci *.combined.gtf *.stats "$bam_out"
    fi
}

# Cufflinks function for multiple SRAs

cufflinks_SRA_multi()
{
    echo "###################"
    echo "Running Cufflinks"
    echo "####################"
    echo "sambamba-0.6.8-linux-static sort --tmpdir=temp -t $num_threads -o "$bam_out"/$f.sorted.bam $f.bam"
    sambamba-0.6.8-linux-static sort --tmpdir=temp -t $num_threads -o "$bam_out"/$f.sorted.bam $f.bam
    rm -r temp $f.bam
    echo "cufflinks "$bam_out"/$f.sorted.bam -p $num_threads -g $referenceannotation -o "$bam_out""
    cufflinks "$bam_out"/$f.sorted.bam -p $num_threads -g $referenceannotation -o "$bam_out"
    mv "$bam_out"/transcripts.gtf "$bam_out"/$f.gtf
    mv "$bam_out"/isoforms.fpkm_tracking "$bam_out"/$f.isoforms.fpkm_tracking
    mv "$bam_out"/skipped.gtf  "$bam_out"/$f.skipped.gtf
    coverge_cuffoff_SRA_multi
    fpkm_cuffoff_SRA_multi
    var=$(grep -vFf listtoremove.txt "$bam_out"/$f.gtf.filtered.gtf | wc -l)
    if [ "$var" -eq 0 ]; then
      rm listtoremove.txt "$bam_out"/$f.gtf.filtered.gtf
    else
      echo "###################"
      echo "Running Cuffcompare"
      echo "####################"
      echo "cuffcompare "$bam_out"/$f.gtf.filtered.gtf -r $referenceannotation -o $f"
      cuffcompare "$bam_out"/$f.gtf.filtered.gtf -r $referenceannotation -o $f
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
        echo "###################"
        echo "Running Cuffmerge"
        echo "####################"
        echo "cuffmerge -o "$bam_out"/merged_out -g $referenceannotation "$bam_out"/gtf_file.txt -p $num_threads"
        cuffmerge -o "$bam_out"/merged_out -g $referenceannotation "$bam_out"/gtf_file.txt -p $num_threads
        rm "$bam_out"/gtf_file.txt
        mkdir index && mv $fbname* index
        mv mapped.txt metrics.txt "$bam_out"
        echo "##############################"
        echo "Pipeline executed successfully"
        echo "##############################"
      else
        echo "Cuffmerge needs more than 1 filtered.gtf file!!!" 1>&2
      fi
    else
      echo "No filtered.gtf files are found for cuffmerge!!!" 1>&2
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
  echo "####################################"
  echo "Starting Reference genome build"
  echo "####################################"
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
            echo "------------------------------------------" >> mapped.txt
            echo "### Mapping percentages of" $filename3 "###" >> mapped.txt
            echo "------------------------------------------" >> mapped.txt
            echo "#######################"
            echo "Running Hisat2 mapping"
            echo "#######################"
            echo "hisat2 -x $fbname --rna-strandness $lib_type -1 ${filename}.fq.gz -2 ${filename2}.fq.gz -p $num_threads -5 $five_trim -3 $three_trim --min-intronlen $min_intl --max-intronlen $max_intl --dta --phred33 -S $filename3.sam --met-file metrics.txt 2>> mapped.txt"
            hisat2 -x $fbname --rna-strandness $lib_type -1 ${filename}.fq.gz -2 ${filename2}.fq.gz -p $num_threads -5 $five_trim -3 $three_trim --min-intronlen $min_intl --max-intronlen $max_intl --dta --phred33 -S $filename3.sam --met-file metrics.txt 2>> mapped.txt 
            echo "samtools view -Sbo $filename3.bam $filename3.sam"
            echo "samtools view -Sbo $filename3.bam $filename3.sam"
            samtools view -Sbo $filename3.bam $filename3.sam
            rm -r $filename3.sam
            mv mapped.txt metrics.txt "$bam_out"
            echo "##################"
            echo "Running Stringtie"
            echo "##################"
            stringtie_non_SRA
            mkdir index && mv $fbname* index
            echo "##############################"
            echo "Pipeline executed successfully"
            echo "##############################"
            exit
        elif [ $numb -eq 1 ] && [ "$cuff_merge" != 0 ]; then
            echo "------------------------------------------" >> mapped.txt
            echo "### Mapping percentages of" $filename3 "###" >> mapped.txt
            echo "------------------------------------------" >> mapped.txt
            echo "#######################"
            echo "Running Hisat2 mapping"
            echo "#######################"
            echo "hisat2 -x $fbname --rna-strandness $lib_type -1 ${filename}.fq.gz -2 ${filename2}.fq.gz -p $num_threads -5 $five_trim -3 $three_trim --min-intronlen $min_intl --max-intronlen $max_intl --dta --phred33 -S $filename3.sam --met-file metrics.txt 2>> mapped.txt"
            hisat2 -x $fbname --rna-strandness $lib_type -1 ${filename}.fq.gz -2 ${filename2}.fq.gz -p $num_threads -5 $five_trim -3 $three_trim --min-intronlen $min_intl --max-intronlen $max_intl --dta --phred33 -S $filename3.sam --met-file metrics.txt 2>> mapped.txt
            echo "samtools view -Sbo $filename3.bam $filename3.sam"
            samtools view -Sbo $filename3.bam $filename3.sam
            rm $filename3.sam
            mv mapped.txt metrics.txt "$bam_out"
            stringtie_non_SRA
            echo "Cuffmerge only works with more than one file!!!"
            mkdir index && mv $fbname* index
            echo "##############################"
            echo "Pipeline executed successfully"
            echo "##############################"
            exit
        elif [ $numb -gt 1 ] && [ "$cuff_merge" == 0 ]; then
            echo "------------------------------------------" >> mapped.txt
            echo "### Mapping percentages of" $filename3 "###">> mapped.txt
            echo "------------------------------------------" >> mapped.txt
            echo "#######################"
            echo "Running Hisat2 mapping"
            echo "#######################"
            echo "hisat2 -x $fbname --rna-strandness $lib_type -1 ${filename}.fq.gz -2 ${filename2}.fq.gz -p $num_threads -5 $five_trim -3 $three_trim --min-intronlen $min_intl --max-intronlen $max_intl --dta --phred33 -S $filename3.sam --met-file metrics.txt 2>> mapped.txt"
            hisat2 -x $fbname --rna-strandness $lib_type -1 ${filename}.fq.gz -2 ${filename2}.fq.gz -p $num_threads -5 $five_trim -3 $three_trim --min-intronlen $min_intl --max-intronlen $max_intl --dta --phred33 -S $filename3.sam --met-file metrics.txt 2>> mapped.txt
            echo "samtools view -Sbo $filename3.bam $filename3.sam"
            samtools view -Sbo $filename3.bam $filename3.sam
            rm $filename3.sam 
            stringtie_non_SRA             
        elif [ $numb -gt 1 ] && [ "$cuff_merge" != 0 ]; then
            echo "------------------------------------------" >> mapped.txt
            echo "### Mapping percentages of" $filename3 "###" >> mapped.txt
            echo "------------------------------------------" >> mapped.txt
            echo "#######################"
            echo "Running Hisat2 mapping"
            echo "#######################"
            echo "hisat2 -x $fbname --rna-strandness $lib_type -1 ${filename}.fq.gz -2 ${filename2}.fq.gz -p $num_threads -5 $five_trim -3 $three_trim --min-intronlen $min_intl --max-intronlen $max_intl --dta --phred33 -S $filename3.sam --met-file metrics.txt 2>> mapped.txt"
            hisat2 -x $fbname --rna-strandness $lib_type -1 ${filename}.fq.gz -2 ${filename2}.fq.gz -p $num_threads -5 $five_trim -3 $three_trim --min-intronlen $min_intl --max-intronlen $max_intl --dta --phred33 -S $filename3.sam --met-file metrics.txt 2>> mapped.txt
            echo "samtools view -Sbo $filename3.bam $filename3.sam"
            samtools view -Sbo $filename3.bam $filename3.sam
            rm $filename3.sam
            stringtie_non_SRA
        fi

      elif [[ "$extension" =~ "fastq.gz" ]]; then
        filename=$(basename "$f" ".fastq.gz")
        filename2=${filename/_R1/_R2}
	      filename3=$(echo $filename | sed 's/.r1//')
        if [ $numb -eq 1 ] && [ "$cuff_merge" == 0 ]; then
            echo "------------------------------------------" >> mapped.txt
            echo "### Mapping percentages of" $filename3 "###" >> mapped.txt
            echo "------------------------------------------" >> mapped.txt
            echo "#######################"
            echo "Running Hisat2 mapping"
            echo "#######################"
            echo "hisat2 -x $fbname --rna-strandness $lib_type -1 ${filename}.fastq.gz -2 ${filename2}.fastq.gz -p $num_threads -5 $five_trim -3 $three_trim --min-intronlen $min_intl --max-intronlen $max_intl --dta --phred33 -S $filename3.sam --met-file metrics.txt 2>> mapped.txt"
            hisat2 -x $fbname --rna-strandness $lib_type -1 ${filename}.fastq.gz -2 ${filename2}.fastq.gz -p $num_threads -5 $five_trim -3 $three_trim --min-intronlen $min_intl --max-intronlen $max_intl --dta --phred33 -S $filename3.sam --met-file metrics.txt 2>> mapped.txt
            echo "samtools view -Sbo $filename3.bam $filename3.sam"
            samtools view -Sbo $filename3.bam $filename3.sam
            rm $filename3.sam
            mv mapped.txt metrics.txt "$bam_out"
            stringtie_non_SRA
            mkdir index && mv $fbname* index
            echo "##############################"
            echo "Pipeline executed successfully"
            echo "##############################"
            exit
        elif [ $numb -eq 1 ] && [ "$cuff_merge" != 0 ]; then
            echo "------------------------------------------" >> mapped.txt
            echo "### Mapping percentages of" $filename3 "###" >> mapped.txt
            echo "------------------------------------------" >> mapped.txt
            echo "#######################"
            echo "Running Hisat2 mapping"
            echo "#######################"
            echo "hisat2 -x $fbname --rna-strandness $lib_type -1 ${filename}.fastq.gz -2 ${filename2}.fastq.gz -p $num_threads -5 $five_trim -3 $three_trim --min-intronlen $min_intl --max-intronlen $max_intl --dta --phred33 -S $filename3.sam --met-file metrics.txt 2>> mapped.txt"
            hisat2 -x $fbname --rna-strandness $lib_type -1 ${filename}.fastq.gz -2 ${filename2}.fastq.gz -p $num_threads -5 $five_trim -3 $three_trim --min-intronlen $min_intl --max-intronlen $max_intl --dta --phred33 -S $filename3.sam --met-file metrics.txt 2>> mapped.txt
            echo "samtools view -Sbo $filename3.bam $filename3.sam"
            samtools view -Sbo $filename3.bam $filename3.sam
            rm $filename3.sam
            mv mapped.txt metrics.txt "$bam_out"
            stringtie_non_SRA
            echo "Cuffmerge only works with more than one file!!!"
            mkdir index && mv $fbname* index
            echo "##############################"
            echo "Pipeline executed successfully"
            echo "##############################"
            exit    
        elif [ $numb -gt 1 ] && [ "$cuff_merge" == 0 ]; then
            echo "------------------------------------------" >> mapped.txt
            echo "### Mapping percentages of" $filename3 "###" >> mapped.txt
            echo "------------------------------------------" >> mapped.txt
            echo "#######################"
            echo "Running Hisat2 mapping"
            echo "#######################"
            echo "hisat2 -x $fbname --rna-strandness $lib_type -1 ${filename}.fastq.gz -2 ${filename2}.fastq.gz -p $num_threads -5 $five_trim -3 $three_trim --min-intronlen $min_intl --max-intronlen $max_intl --dta --phred33 -S $filename3.sam --met-file metrics.txt 2>> mapped.txt"
            hisat2 -x $fbname --rna-strandness $lib_type -1 ${filename}.fastq.gz -2 ${filename2}.fastq.gz -p $num_threads -5 $five_trim -3 $three_trim --min-intronlen $min_intl --max-intronlen $max_intl --dta --phred33 -S $filename3.sam --met-file metrics.txt 2>> mapped.txt
            echo "samtools view -Sbo $filename3.bam $filename3.sam"
            samtools view -Sbo $filename3.bam $filename3.sam
            rm $filename3.sam
            stringtie_non_SRA
        elif [ $numb -gt 1 ] && [ "$cuff_merge" != 0 ]; then
            echo "------------------------------------------" >> mapped.txt
            echo "### Mapping percentages of" $filename3 "###" >> mapped.txt
            echo "------------------------------------------" >> mapped.txt
            echo "#######################"
            echo "Running Hisat2 mapping"
            echo "#######################"
            echo "hisat2 -x $fbname --rna-strandness $lib_type -1 ${filename}.fastq.gz -2 ${filename2}.fastq.gz -p $num_threads -5 $five_trim -3 $three_trim --min-intronlen $min_intl --max-intronlen $max_intl --dta --phred33 -S $filename3.sam --met-file metrics.txt 2>> mapped.txt"
            hisat2 -x $fbname --rna-strandness $lib_type -1 ${filename}.fastq.gz -2 ${filename2}.fastq.gz -p $num_threads -5 $five_trim -3 $three_trim --min-intronlen $min_intl --max-intronlen $max_intl --dta --phred33 -S $filename3.sam --met-file metrics.txt 2>> mapped.txt
            echo "samtools view -Sbo $filename3.bam $filename3.sam"
            samtools view -Sbo $filename3.bam $filename3.sam
            rm $filename3.sam
            stringtie_non_SRA
        fi

      elif [[ "$extension" =~ "fq" ]]; then
        filename=$(basename "$f" ".fq")
        filename2=${filename/_R1/_R2}
	      filename3=$(echo $filename | sed 's/_R1//')
        if [ $numb -eq 1 ] && [ "$cuff_merge" == 0 ]; then
            echo "------------------------------------------" >> mapped.txt
            echo "### Mapping percentages of" $filename3 "###" >> mapped.txt
            echo "------------------------------------------" >> mapped.txt
            echo "#######################"
            echo "Running Hisat2 mapping"
            echo "#######################"
            echo "hisat2 -x $fbname --rna-strandness $lib_type -1 ${filename}.fq -2 ${filename2}.fq -p $num_threads -5 $five_trim -3 $three_trim --min-intronlen $min_intl --max-intronlen $max_intl --dta --phred33 -S $filename3.sam --met-file metrics.txt 2>> mapped.txt"
            hisat2 -x $fbname --rna-strandness $lib_type -1 ${filename}.fq -2 ${filename2}.fq -p $num_threads -5 $five_trim -3 $three_trim --min-intronlen $min_intl --max-intronlen $max_intl --dta --phred33 -S $filename3.sam --met-file metrics.txt 2>> mapped.txt
            echo "samtools view -Sbo $filename3.bam $filename3.sam"
            samtools view -Sbo $filename3.bam $filename3.sam
            rm $filename3.sam
            mv mapped.txt metrics.txt "$bam_out"
            stringtie_non_SRA
            mkdir index && mv $fbname* index
            echo "##############################"
            echo "Pipeline executed successfully"
            echo "##############################"
            exit
        elif [ $numb -eq 1 ] && [ "$cuff_merge" != 0 ]; then
            echo "------------------------------------------" >> mapped.txt
            echo "### Mapping percentages of" $filename3 "###" >> mapped.txt
            echo "------------------------------------------" >> mapped.txt
            echo "#######################"
            echo "Running Hisat2 mapping"
            echo "#######################"
            echo "hisat2 -x $fbname --rna-strandness $lib_type -1 ${filename}.fq -2 ${filename2}.fq -p $num_threads -5 $five_trim -3 $three_trim --min-intronlen $min_intl --max-intronlen $max_intl --dta --phred33 -S $filename3.sam --met-file metrics.txt 2>> mapped.txt"
            hisat2 -x $fbname --rna-strandness $lib_type -1 ${filename}.fq -2 ${filename2}.fq -p $num_threads -5 $five_trim -3 $three_trim --min-intronlen $min_intl --max-intronlen $max_intl --dta --phred33 -S $filename3.sam --met-file metrics.txt 2>> mapped.txt
            echo "samtools view -Sbo $filename3.bam $filename3.sam"
            samtools view -Sbo $filename3.bam $filename3.sam
            rm $filename3.sam
            mv mapped.txt metrics.txt "$bam_out"
            stringtie_non_SRA
            mkdir index && mv $fbname* index
            echo "Cuffmerge only works with more than one file!!!"
            echo "##############################"
            echo "Pipeline executed successfully"
            echo "##############################"
            exit
        elif [ $numb -gt 1 ] && [ "$cuff_merge" == 0 ]; then
            echo "------------------------------------------" >> mapped.txt
            echo "### Mapping percentages of" $filename3 "###" >> mapped.txt
            echo "------------------------------------------" >> mapped.txt
            echo "#######################"
            echo "Running Hisat2 mapping"
            echo "#######################"
            echo "hisat2 -x $fbname --rna-strandness $lib_type -1 ${filename}.fq -2 ${filename2}.fq -p $num_threads -5 $five_trim -3 $three_trim --min-intronlen $min_intl --max-intronlen $max_intl --dta --phred33 -S $filename3.sam --met-file metrics.txt 2>> mapped.txt"
            hisat2 -x $fbname --rna-strandness $lib_type -1 ${filename}.fq -2 ${filename2}.fq -p $num_threads -5 $five_trim -3 $three_trim --min-intronlen $min_intl --max-intronlen $max_intl --dta --phred33 -S $filename3.sam --met-file metrics.txt 2>> mapped.txt
            echo "samtools view -Sbo $filename3.bam $filename3.sam"
            samtools view -Sbo $filename3.bam $filename3.sam
            rm $filename3.sam
            stringtie_non_SRA   
        elif [ $numb -gt 1 ] && [ "$cuff_merge" != 0 ]; then
            echo "------------------------------------------" >> mapped.txt
            echo "### Mapping percentages of" $filename3 "###" >> mapped.txt
            echo "------------------------------------------" >> mapped.txt
            echo "#######################"
            echo "Running Hisat2 mapping"
            echo "#######################"
            echo "hisat2 -x $fbname --rna-strandness $lib_type -1 ${filename}.fq -2 ${filename2}.fq -p $num_threads -5 $five_trim -3 $three_trim --min-intronlen $min_intl --max-intronlen $max_intl --dta --phred33 -S $filename3.sam --met-file metrics.txt 2>> mapped.txt"
            hisat2 -x $fbname --rna-strandness $lib_type -1 ${filename}.fq -2 ${filename2}.fq -p $num_threads -5 $five_trim -3 $three_trim --min-intronlen $min_intl --max-intronlen $max_intl --dta --phred33 -S $filename3.sam --met-file metrics.txt 2>> mapped.txt
            echo "samtools view -Sbo $filename3.bam $filename3.sam"
            samtools view -Sbo $filename3.bam $filename3.sam
            rm $filename3.sam
            stringtie_non_SRA
        fi

      elif [[ "$extension" =~ "fastq" ]]; then
        filename=$(basename "$f" ".fastq")
        filename2=${filename/_R1/_R2}
	      filename3=$(echo $filename | sed 's/_R1//')
        if [ $numb -eq 1 ] && [ "$cuff_merge" == 0 ]; then
            echo "------------------------------------------" >> mapped.txt
            echo "### Mapping percentages of" $filename3 "###" >> mapped.txt
            echo "------------------------------------------" >> mapped.txt
            echo "#######################"
            echo "Running Hisat2 mapping"
            echo "#######################"
            echo "hisat2 -x $fbname --rna-strandness $lib_type -1 ${filename}.fastq -2 ${filename2}.fastq -p $num_threads -5 $five_trim -3 $three_trim --min-intronlen $min_intl --max-intronlen $max_intl --dta --phred33 -S $filename3.sam --met-file metrics.txt 2>> mapped.txt"
            hisat2 -x $fbname --rna-strandness $lib_type -1 ${filename}.fastq -2 ${filename2}.fastq -p $num_threads -5 $five_trim -3 $three_trim --min-intronlen $min_intl --max-intronlen $max_intl --dta --phred33 -S $filename3.sam --met-file metrics.txt 2>> mapped.txt
            echo "samtools view -Sbo $filename3.bam $filename3.sam"
            samtools view -Sbo $filename3.bam $filename3.sam
            rm $filename3.sam
            mv mapped.txt metrics.txt "$bam_out"
            stringtie_non_SRA
            mkdir index && mv $fbname* index
            echo "##############################"
            echo "Pipeline executed successfully"
            echo "##############################"
            exit
        elif [ $numb -eq 1 ] && [ "$cuff_merge" != 0 ]; then
            echo "------------------------------------------" >> mapped.txt
            echo "### Mapping percentages of" $filename3 "###" >> mapped.txt
            echo "------------------------------------------" >> mapped.txt
            echo "#######################"
            echo "Running Hisat2 mapping"
            echo "#######################"
            echo "hisat2 -x $fbname --rna-strandness $lib_type -1 ${filename}.fastq -2 ${filename2}.fastq -p $num_threads -5 $five_trim -3 $three_trim --min-intronlen $min_intl --max-intronlen $max_intl --dta --phred33 -S $filename3.sam --met-file metrics.txt 2>> mapped.txt"
            hisat2 -x $fbname --rna-strandness $lib_type -1 ${filename}.fastq -2 ${filename2}.fastq -p $num_threads -5 $five_trim -3 $three_trim --min-intronlen $min_intl --max-intronlen $max_intl --dta --phred33 -S $filename3.sam --met-file metrics.txt 2>> mapped.txt
            echo "samtools view -Sbo $filename3.bam $filename3.sam"
            samtools view -Sbo $filename3.bam $filename3.sam
            rm $filename3.sam
            mv mapped.txt metrics.txt "$bam_out"
            stringtie_non_SRA
            mkdir index && mv $fbname* index
            echo "Cuffmerge only works with more than one file!!!"
            echo "##############################"
            echo "Pipeline executed successfully"
            echo "##############################"
            exit
        elif [ $numb -gt 1 ] && [ "$cuff_merge" == 0 ]; then
            echo "------------------------------------------" >> mapped.txt
            echo "### Mapping percentages of" $filename3 "###" >> mapped.txt
            echo "------------------------------------------" >> mapped.txt
            echo "#######################"
            echo "Running Hisat2 mapping"
            echo "#######################"
            echo "hisat2 -x $fbname --rna-strandness $lib_type -1 ${filename}.fastq -2 ${filename2}.fastq -p $num_threads -5 $five_trim -3 $three_trim --min-intronlen $min_intl --max-intronlen $max_intl --dta --phred33 -S $filename3.sam --met-file metrics.txt 2>> mapped.txt"
            hisat2 -x $fbname --rna-strandness $lib_type -1 ${filename}.fastq -2 ${filename2}.fastq -p $num_threads -5 $five_trim -3 $three_trim --min-intronlen $min_intl --max-intronlen $max_intl --dta --phred33 -S $filename3.sam --met-file metrics.txt 2>> mapped.txt
            echo "samtools view -Sbo $filename3.bam $filename3.sam"
            samtools view -Sbo $filename3.bam $filename3.sam
            rm $filename3.sam
            stringtie_non_SRA               
        elif [ $numb -gt 1 ] && [ "$cuff_merge" != 0 ]; then
            echo "------------------------------------------" >> mapped.txt
            echo "### Mapping percentages of" $filename3 "###" >> mapped.txt
            echo "------------------------------------------" >> mapped.txt
            echo "#######################"
            echo "Running Hisat2 mapping"
            echo "#######################"
            echo "hisat2 -x $fbname --rna-strandness $lib_type -1 ${filename}.fastq -2 ${filename2}.fastq -p $num_threads -5 $five_trim -3 $three_trim --min-intronlen $min_intl --max-intronlen $max_intl --dta --phred33 -S $filename3.sam --met-file metrics.txt 2>> mapped.txt"
            hisat2 -x $fbname --rna-strandness $lib_type -1 ${filename}.fastq -2 ${filename2}.fastq -p $num_threads -5 $five_trim -3 $three_trim --min-intronlen $min_intl --max-intronlen $max_intl --dta --phred33 -S $filename3.sam --met-file metrics.txt 2>> mapped.txt
            echo "samtools view -Sbo $filename3.bam $filename3.sam"
            samtools view -Sbo $filename3.bam $filename3.sam
            rm $filename3.sam
            stringtie_non_SRA
        fi   

      elif [ "$extension" != "fastq" ] || [ "$extension" != "fq" ] || [ "$extension" != "fastq.gz" ] || [ "$extension" != "fq.gz" ]; then
        echo "The extension" "$extension" "is not supported. Only .fq, .fq.gz, .fastq, .fastq.gz are only supported" 1>&2        
        exit 64
      fi 
    done
    
    if [ $numb -gt 1 ] && [ "$cuff_merge" == 0 ]; then
      mkdir index && mv $fbname* index
      mv mapped.txt metrics.txt "$bam_out"
      echo "##############################"
      echo "Pipeline executed successfully"
      echo "##############################"
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
            echo "------------------------------------------" >> mapped.txt
            echo "### Mapping percentages of" $filename3 "###" >> mapped.txt
            echo "------------------------------------------" >> mapped.txt
            echo "#######################"
            echo "Running Hisat2 mapping"
            echo "#######################"
            echo "hisat2 -x $fbname --rna-strandness $lib_type -1 ${filename}.fq.gz -2 ${filename2}.fq.gz -p $num_threads -5 $five_trim -3 $three_trim --min-intronlen $min_intl --max-intronlen $max_intl --dta --phred64 -S $filename3.sam --met-file metrics.txt 2>> mapped.txt"
            hisat2 -x $fbname --rna-strandness $lib_type -1 ${filename}.fq.gz -2 ${filename2}.fq.gz -p $num_threads -5 $five_trim -3 $three_trim --min-intronlen $min_intl --max-intronlen $max_intl --dta --phred64 -S $filename3.sam --met-file metrics.txt 2>> mapped.txt
            echo "samtools view -Sbo $filename3.bam $filename3.sam"
            samtools view -Sbo $filename3.bam $filename3.sam
            rm $filename3.sam
            mv mapped.txt metrics.txt "$bam_out"
            stringtie_non_SRA
            mkdir index && mv $fbname* index
            echo "##############################"
            echo "Pipeline executed successfully"
            echo "##############################"
            exit
        elif [ $numb -eq 1 ] && [ "$cuff_merge" != 0 ]; then
            echo "------------------------------------------" >> mapped.txt
            echo "### Mapping percentages of" $filename3 "###" >> mapped.txt
            echo "------------------------------------------" >> mapped.txt
            echo "#######################"
            echo "Running Hisat2 mapping"
            echo "#######################"
            echo "hisat2 -x $fbname --rna-strandness $lib_type -1 ${filename}.fq.gz -2 ${filename2}.fq.gz -p $num_threads -5 $five_trim -3 $three_trim --min-intronlen $min_intl --max-intronlen $max_intl --dta --phred64 -S $filename3.sam --met-file metrics.txt 2>> mapped.txt"
            hisat2 -x $fbname --rna-strandness $lib_type -1 ${filename}.fq.gz -2 ${filename2}.fq.gz -p $num_threads -5 $five_trim -3 $three_trim --min-intronlen $min_intl --max-intronlen $max_intl --dta --phred64 -S $filename3.sam --met-file metrics.txt 2>> mapped.txt
            echo "samtools view -Sbo $filename3.bam $filename3.sam"
            samtools view -Sbo $filename3.bam $filename3.sam
            rm $filename3.sam
            mv mapped.txt metrics.txt "$bam_out"
            stringtie_non_SRA
            mkdir index && mv $fbname* index
            echo "Cuffmerge only works with more than one file!!!"
            echo "##############################"
            echo "Pipeline executed successfully"
            echo "##############################"
            exit    
        elif [ $numb -gt 1 ] && [ "$cuff_merge" == 0 ]; then
            echo "------------------------------------------" >> mapped.txt
            echo "### Mapping percentages of" $filename3 "###" >> mapped.txt
            echo "------------------------------------------" >> mapped.txt
            echo "#######################"
            echo "Running Hisat2 mapping"
            echo "#######################"
            echo "hisat2 -x $fbname --rna-strandness $lib_type -1 ${filename}.fq.gz -2 ${filename2}.fq.gz -p $num_threads -5 $five_trim -3 $three_trim --min-intronlen $min_intl --max-intronlen $max_intl --dta --phred64 -S $filename3.sam --met-file metrics.txt 2>> mapped.txt"
            hisat2 -x $fbname --rna-strandness $lib_type -1 ${filename}.fq.gz -2 ${filename2}.fq.gz -p $num_threads -5 $five_trim -3 $three_trim --min-intronlen $min_intl --max-intronlen $max_intl --dta --phred64 -S $filename3.sam --met-file metrics.txt 2>> mapped.txt
            echo "samtools view -Sbo $filename3.bam $filename3.sam"
            samtools view -Sbo $filename3.bam $filename3.sam
            rm $filename3.sam
            stringtie_non_SRA
        elif [ $numb -gt 1 ] && [ "$cuff_merge" != 0 ]; then
            echo "------------------------------------------" >> mapped.txt
            echo "### Mapping percentages of" $filename3 "###" >> mapped.txt
            echo "------------------------------------------" >> mapped.txt
            echo "#######################"
            echo "Running Hisat2 mapping"
            echo "#######################"
            echo "hisat2 -x $fbname --rna-strandness $lib_type -1 ${filename}.fq.gz -2 ${filename2}.fq.gz -p $num_threads -5 $five_trim -3 $three_trim --min-intronlen $min_intl --max-intronlen $max_intl --dta --phred64 -S $filename3.sam --met-file metrics.txt 2>> mapped.txt"
            hisat2 -x $fbname --rna-strandness $lib_type -1 ${filename}.fq.gz -2 ${filename2}.fq.gz -p $num_threads -5 $five_trim -3 $three_trim --min-intronlen $min_intl --max-intronlen $max_intl --dta --phred64 -S $filename3.sam --met-file metrics.txt 2>> mapped.txt
            echo "samtools view -Sbo $filename3.bam $filename3.sam"
            samtools view -Sbo $filename3.bam $filename3.sam
            rm $filename3.sam
            stringtie_non_SRA
        fi

      elif [[ "$extension" =~ "fastq.gz" ]]; then
        filename=$(basename "$f" ".fastq.gz")
        filename2=${filename/_R1/_R2}
	      filename3=$(echo $filename | sed 's/_R1//')
        if [ $numb -eq 1 ] && [ "$cuff_merge" == 0 ]; then
            echo "------------------------------------------" >> mapped.txt
            echo "### Mapping percentages of" $filename3 "###" >> mapped.txt
            echo "------------------------------------------" >> mapped.txt
            echo "#######################"
            echo "Running Hisat2 mapping"
            echo "#######################"
            echo "hisat2 -x $fbname --rna-strandness $lib_type -1 ${filename}.fastq.gz -2 ${filename2}.fastq.gz -p $num_threads -5 $five_trim -3 $three_trim --min-intronlen $min_intl --max-intronlen $max_intl --dta --phred64 -S $filename3.sam --met-file metrics.txt 2>> mapped.txt"
            hisat2 -x $fbname --rna-strandness $lib_type -1 ${filename}.fastq.gz -2 ${filename2}.fastq.gz -p $num_threads -5 $five_trim -3 $three_trim --min-intronlen $min_intl --max-intronlen $max_intl --dta --phred64 -S $filename3.sam --met-file metrics.txt 2>> mapped.txt
            echo "samtools view -Sbo $filename3.bam $filename3.sam"
            samtools view -Sbo $filename3.bam $filename3.sam
            rm $filename3.sam
            mv mapped.txt metrics.txt "$bam_out"
            stringtie_non_SRA
            mkdir index && mv $fbname* index
            echo "##############################"
            echo "Pipeline executed successfully"
            echo "##############################"
            exit
        elif [ $numb -eq 1 ] && [ "$cuff_merge" != 0 ]; then
            echo "------------------------------------------" >> mapped.txt
            echo "### Mapping percentages of" $filename3 "###" >> mapped.txt
            echo "------------------------------------------" >> mapped.txt
            echo "#######################"
            echo "Running Hisat2 mapping"
            echo "#######################"
            echo "hisat2 -x $fbname --rna-strandness $lib_type -1 ${filename}.fastq.gz -2 ${filename2}.fastq.gz -p $num_threads -5 $five_trim -3 $three_trim --min-intronlen $min_intl --max-intronlen $max_intl --dta --phred64 -S $filename3.sam --met-file metrics.txt 2>> mapped.txt"
            hisat2 -x $fbname --rna-strandness $lib_type -1 ${filename}.fastq.gz -2 ${filename2}.fastq.gz -p $num_threads -5 $five_trim -3 $three_trim --min-intronlen $min_intl --max-intronlen $max_intl --dta --phred64 -S $filename3.sam --met-file metrics.txt 2>> mapped.txt
            echo "samtools view -Sbo $filename3.bam $filename3.sam"
            samtools view -Sbo $filename3.bam $filename3.sam
            rm $filename3.sam
            mv mapped.txt metrics.txt "$bam_out"
            stringtie_non_SRA
            mkdir index && mv $fbname* index
            echo "Cuffmerge only works with more than one file!!!"
            echo "##############################"
            echo "Pipeline executed successfully"
            echo "##############################"
            exit    
        elif [ $numb -gt 1 ] && [ "$cuff_merge" == 0 ]; then
            echo "------------------------------------------" >> mapped.txt
            echo "### Mapping percentages of" $filename3 "###" >> mapped.txt
            echo "------------------------------------------" >> mapped.txt
            echo "#######################"
            echo "Running Hisat2 mapping"
            echo "#######################"
            echo "hisat2 -x $fbname --rna-strandness $lib_type -1 ${filename}.fastq.gz -2 ${filename2}.fastq.gz -p $num_threads -5 $five_trim -3 $three_trim --min-intronlen $min_intl --max-intronlen $max_intl --dta --phred64 -S $filename3.sam --met-file metrics.txt 2>> mapped.txt"
            hisat2 -x $fbname --rna-strandness $lib_type -1 ${filename}.fastq.gz -2 ${filename2}.fastq.gz -p $num_threads -5 $five_trim -3 $three_trim --min-intronlen $min_intl --max-intronlen $max_intl --dta --phred64 -S $filename3.sam --met-file metrics.txt 2>> mapped.txt
            echo "samtools view -Sbo $filename3.bam $filename3.sam"
            samtools view -Sbo $filename3.bam $filename3.sam
            rm $filename3.sam
            stringtie_non_SRA
        elif [ $numb -gt 1 ] && [ "$cuff_merge" != 0 ]; then
            echo "------------------------------------------" >> mapped.txt
            echo "### Mapping percentages of" $filename3 "###" >> mapped.txt
            echo "------------------------------------------" >> mapped.txt
            echo "#######################"
            echo "Running Hisat2 mapping"
            echo "#######################"
            echo "hisat2 -x $fbname --rna-strandness $lib_type -1 ${filename}.fastq.gz -2 ${filename2}.fastq.gz -p $num_threads -5 $five_trim -3 $three_trim --min-intronlen $min_intl --max-intronlen $max_intl --dta --phred64 -S $filename3.sam --met-file metrics.txt 2>> mapped.txt"
            hisat2 -x $fbname --rna-strandness $lib_type -1 ${filename}.fastq.gz -2 ${filename2}.fastq.gz -p $num_threads -5 $five_trim -3 $three_trim --min-intronlen $min_intl --max-intronlen $max_intl --dta --phred64 -S $filename3.sam --met-file metrics.txt 2>> mapped.txt
            echo "samtools view -Sbo $filename3.bam $filename3.sam"
            samtools view -Sbo $filename3.bam $filename3.sam
            rm $filename3.sam
            stringtie_non_SRA
        fi

      elif [[ "$extension" =~ "fq" ]]; then
        filename=$(basename "$f" ".fq")
        filename2=${filename/_R1/_R2}
	      filename3=$(echo $filename | sed 's/_R1//')

        if [ $numb -eq 1 ] && [ "$cuff_merge" == 0 ]; then
            echo "------------------------------------------" >> mapped.txt
            echo "### Mapping percentages of" $filename3 "###" >> mapped.txt
            echo "------------------------------------------" >> mapped.txt
            echo "#######################"
            echo "Running Hisat2 mapping"
            echo "#######################"
            echo "hisat2 -x $fbname --rna-strandness $lib_type -1 ${filename}.fq -2 ${filename2}.fq -p $num_threads -5 $five_trim -3 $three_trim --min-intronlen $min_intl --max-intronlen $max_intl --dta --phred64 -S $filename3.sam --met-file metrics.txt 2>> mapped.txt"
            hisat2 -x $fbname --rna-strandness $lib_type -1 ${filename}.fq -2 ${filename2}.fq -p $num_threads -5 $five_trim -3 $three_trim --min-intronlen $min_intl --max-intronlen $max_intl --dta --phred64 -S $filename3.sam --met-file metrics.txt 2>> mapped.txt
            echo "samtools view -Sbo $filename3.bam $filename3.sam"
            samtools view -Sbo $filename3.bam $filename3.sam
            rm $filename3.sam
            mv mapped.txt metrics.txt "$bam_out"
            stringtie_non_SRA
            mkdir index && mv $fbname* index
            echo "##############################"
            echo "Pipeline executed successfully"
            echo "##############################"
            exit
        elif [ $numb -eq 1 ] && [ "$cuff_merge" != 0 ]; then
            echo "------------------------------------------" >> mapped.txt
            echo "### Mapping percentages of" $filename3 "###" >> mapped.txt
            echo "------------------------------------------" >> mapped.txt
            echo "#######################"
            echo "Running Hisat2 mapping"
            echo "#######################"
            echo "hisat2 -x $fbname --rna-strandness $lib_type -1 ${filename}.fq -2 ${filename2}.fq -p $num_threads -5 $five_trim -3 $three_trim --min-intronlen $min_intl --max-intronlen $max_intl --dta --phred64 -S $filename3.sam --met-file metrics.txt 2>> mapped.txt"
            hisat2 -x $fbname --rna-strandness $lib_type -1 ${filename}.fq -2 ${filename2}.fq -p $num_threads -5 $five_trim -3 $three_trim --min-intronlen $min_intl --max-intronlen $max_intl --dta --phred64 -S $filename3.sam --met-file metrics.txt 2>> mapped.txt
            echo "samtools view -Sbo $filename3.bam $filename3.sam"
            samtools view -Sbo $filename3.bam $filename3.sam
            rm $filename3.sam
            mv mapped.txt metrics.txt "$bam_out"
            stringtie_non_SRA
            mkdir index && mv $fbname* index
            echo "Cuffmerge only works with more than one file!!!"
            echo "##############################"
            echo "Pipeline executed successfully"
            echo "##############################"
            exit    
        elif [ $numb -gt 1 ] && [ "$cuff_merge" == 0 ]; then
            echo "------------------------------------------" >> mapped.txt
            echo "### Mapping percentages of" $filename3 "###" >> mapped.txt
            echo "------------------------------------------" >> mapped.txt
            echo "#######################"
            echo "Running Hisat2 mapping"
            echo "#######################"
            echo "hisat2 -x $fbname --rna-strandness $lib_type -1 ${filename}.fq -2 ${filename2}.fq -p $num_threads -5 $five_trim -3 $three_trim --min-intronlen $min_intl --max-intronlen $max_intl --dta --phred64 -S $filename3.sam --met-file metrics.txt 2>> mapped.txt"
            hisat2 -x $fbname --rna-strandness $lib_type -1 ${filename}.fq -2 ${filename2}.fq -p $num_threads -5 $five_trim -3 $three_trim --min-intronlen $min_intl --max-intronlen $max_intl --dta --phred64 -S $filename3.sam --met-file metrics.txt 2>> mapped.txt
            echo "samtools view -Sbo $filename3.bam $filename3.sam"
            samtools view -Sbo $filename3.bam $filename3.sam
            rm $filename3.sam
            stringtie_non_SRA
        elif [ $numb -gt 1 ] && [ "$cuff_merge" != 0 ]; then
            echo "------------------------------------------" >> mapped.txt
            echo "### Mapping percentages of" $filename3 "###" >> mapped.txt
            echo "------------------------------------------" >> mapped.txt
            echo "#######################"
            echo "Running Hisat2 mapping"
            echo "#######################"
            echo "hisat2 -x $fbname --rna-strandness $lib_type -1 ${filename}.fq -2 ${filename2}.fq -p $num_threads -5 $five_trim -3 $three_trim --min-intronlen $min_intl --max-intronlen $max_intl --dta --phred64 -S $filename3.sam --met-file metrics.txt 2>> mapped.txt"
            hisat2 -x $fbname --rna-strandness $lib_type -1 ${filename}.fq -2 ${filename2}.fq -p $num_threads -5 $five_trim -3 $three_trim --min-intronlen $min_intl --max-intronlen $max_intl --dta --phred64 -S $filename3.sam --met-file metrics.txt 2>> mapped.txt
            echo "samtools view -Sbo $filename3.bam $filename3.sam"
            samtools view -Sbo $filename3.bam $filename3.sam
            rm $filename3.sam
            stringtie_non_SRA
        fi

      elif [[ "$extension" =~ "fastq" ]]; then
        filename=$(basename "$f" ".fastq")
        filename2=${filename/_R1/_R2}
	      filename3=$(echo $filename | sed 's/_R1//')

        if [ $numb -eq 1 ] && [ "$cuff_merge" == 0 ]; then
            echo "------------------------------------------" >> mapped.txt
            echo "### Mapping percentages of" $filename3 "###" >> mapped.txt
            echo "------------------------------------------" >> mapped.txt
            echo "#######################"
            echo "Running Hisat2 mapping"
            echo "#######################"
            echo "hisat2 -x $fbname --rna-strandness $lib_type -1 ${filename}.fastq -2 ${filename2}.fastq -p $num_threads -5 $five_trim -3 $three_trim --min-intronlen $min_intl --max-intronlen $max_intl --dta --phred64 -S $filename3.sam --met-file metrics.txt 2>> mapped.txt"
            hisat2 -x $fbname --rna-strandness $lib_type -1 ${filename}.fastq -2 ${filename2}.fastq -p $num_threads -5 $five_trim -3 $three_trim --min-intronlen $min_intl --max-intronlen $max_intl --dta --phred64 -S $filename3.sam --met-file metrics.txt 2>> mapped.txt
            echo "samtools view -Sbo $filename3.bam $filename3.sam"
            samtools view -Sbo $filename3.bam $filename3.sam
            rm $filename3.sam
            mv mapped.txt metrics.txt "$bam_out"
            stringtie_non_SRA
            mkdir index && mv $fbname* index
            echo "##############################"
            echo "Pipeline executed successfully"
            echo "##############################"
            exit
        elif [ $numb -eq 1 ] && [ "$cuff_merge" != 0 ]; then
            echo "------------------------------------------" >> mapped.txt
            echo "### Mapping percentages of" $filename3 "###" >> mapped.txt
            echo "------------------------------------------" >> mapped.txt
            echo "#######################"
            echo "Running Hisat2 mapping"
            echo "#######################"
            echo "hisat2 -x $fbname --rna-strandness $lib_type -1 ${filename}.fastq -2 ${filename2}.fastq -p $num_threads -5 $five_trim -3 $three_trim --min-intronlen $min_intl --max-intronlen $max_intl --dta --phred64 -S $filename3.sam --met-file metrics.txt 2>> mapped.txt"
            hisat2 -x $fbname --rna-strandness $lib_type -1 ${filename}.fastq -2 ${filename2}.fastq -p $num_threads -5 $five_trim -3 $three_trim --min-intronlen $min_intl --max-intronlen $max_intl --dta --phred64 -S $filename3.sam --met-file metrics.txt 2>> mapped.txt
            echo "samtools view -Sbo $filename3.bam $filename3.sam"
            samtools view -Sbo $filename3.bam $filename3.sam
            rm $filename3.sam
            mv mapped.txt metrics.txt "$bam_out"
            stringtie_non_SRA
            mkdir index && mv $fbname* index
            echo "Cuffmerge only works with more than one file!!!"
            echo "##############################"
            echo "Pipeline executed successfully"
            echo "##############################"
            exit    
        elif [ $numb -gt 1 ] && [ "$cuff_merge" == 0 ]; then
            echo "------------------------------------------" >> mapped.txt
            echo "### Mapping percentages of" $filename3 "###" >> mapped.txt
            echo "------------------------------------------" >> mapped.txt
            echo "#######################"
            echo "Running Hisat2 mapping"
            echo "#######################"
            echo "hisat2 -x $fbname --rna-strandness $lib_type -1 ${filename}.fastq -2 ${filename2}.fastq -p $num_threads -5 $five_trim -3 $three_trim --min-intronlen $min_intl --max-intronlen $max_intl --dta --phred64 -S $filename3.sam --met-file metrics.txt 2>> mapped.txt"
            hisat2 -x $fbname --rna-strandness $lib_type -1 ${filename}.fastq -2 ${filename2}.fastq -p $num_threads -5 $five_trim -3 $three_trim --min-intronlen $min_intl --max-intronlen $max_intl --dta --phred64 -S $filename3.sam --met-file metrics.txt 2>> mapped.txt
            echo "samtools view -Sbo $filename3.bam $filename3.sam"
            samtools view -Sbo $filename3.bam $filename3.sam
            rm $filename3.sam
            stringtie_non_SRA
        elif [ $numb -gt 1 ] && [ "$cuff_merge" != 0 ]; then
            echo "------------------------------------------" >> mapped.txt
            echo "### Mapping percentages of" $filename3 "###" >> mapped.txt
            echo "------------------------------------------" >> mapped.txt
            echo "#######################"
            echo "Running Hisat2 mapping"
            echo "#######################"
            echo "hisat2 -x $fbname --rna-strandness $lib_type -1 ${filename}.fastq -2 ${filename2}.fastq -p $num_threads -5 $five_trim -3 $three_trim --min-intronlen $min_intl --max-intronlen $max_intl --dta --phred64 -S $filename3.sam --met-file metrics.txt 2>> mapped.txt"
            hisat2 -x $fbname --rna-strandness $lib_type -1 ${filename}.fastq -2 ${filename2}.fastq -p $num_threads -5 $five_trim -3 $three_trim --min-intronlen $min_intl --max-intronlen $max_intl --dta --phred64 -S $filename3.sam --met-file metrics.txt 2>> mapped.txt
            echo "samtools view -Sbo $filename3.bam $filename3.sam"
            samtools view -Sbo $filename3.bam $filename3.sam
            rm $filename3.sam
            stringtie_non_SRA
        fi   

      elif [ "$extension" != "fastq" ] || [ "$extension" != "fq" ] || [ "$extension" != "fastq.gz" ] || [ "$extension" != "fq.gz" ]; then
        echo "The extension" "$extension" "is not supported. Only .fq, .fq.gz, .fastq, .fastq.gz are only supported!!!!" 1>&2        
        exit 64
      fi 
    done
    
    if [ $numb -gt 1 ] && [ "$cuff_merge" == 0 ]; then
      mkdir index && mv $fbname* index
      mv mapped.txt metrics.txt "$bam_out"
      echo "##############################"
      echo "Pipeline executed successfully"
      echo "##############################"
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
            echo "------------------------------------------" >> mapped.txt
            echo "### Mapping percentages of" $filename3 "###" >> mapped.txt
            echo "------------------------------------------" >> mapped.txt
            echo "#######################"
            echo "Running Hisat2 mapping"
            echo "#######################"
            echo "hisat2 -x $fbname --rna-strandness $lib_type -1 ${filename}.fq.gz -2 ${filename2}.fq.gz -p $num_threads -5 $five_trim -3 $three_trim --min-intronlen $min_intl --max-intronlen $max_intl --dta --phred33 -S $filename3.sam --met-file metrics.txt 2>> mapped.txt"
            hisat2 -x $fbname --rna-strandness $lib_type -1 ${filename}.fq.gz -2 ${filename2}.fq.gz -p $num_threads -5 $five_trim -3 $three_trim --min-intronlen $min_intl --max-intronlen $max_intl --dta --phred33 -S $filename3.sam --met-file metrics.txt 2>> mapped.txt
            echo "samtools view -Sbo $filename3.bam $filename3.sam"
            samtools view -Sbo $filename3.bam $filename3.sam
            rm $filename3.sam
            mv mapped.txt metrics.txt "$bam_out"
            cufflinks_non_SRA
            mkdir index && mv $fbname* index
            echo "##############################"
            echo "Pipeline executed successfully"
            echo "##############################"
            exit
        elif [ $numb -eq 1 ] && [ "$cuff_merge" != 0 ]; then
            echo "------------------------------------------" >> mapped.txt
            echo "### Mapping percentages of" $filename3 "###" >> mapped.txt
            echo "------------------------------------------" >> mapped.txt
            echo "#######################"
            echo "Running Hisat2 mapping"
            echo "#######################"
            echo "hisat2 -x $fbname --rna-strandness $lib_type -1 ${filename}.fq.gz -2 ${filename2}.fq.gz -p $num_threads -5 $five_trim -3 $three_trim --min-intronlen $min_intl --max-intronlen $max_intl --dta --phred33 -S $filename3.sam --met-file metrics.txt 2>> mapped.txt"
            hisat2 -x $fbname --rna-strandness $lib_type -1 ${filename}.fq.gz -2 ${filename2}.fq.gz -p $num_threads -5 $five_trim -3 $three_trim --min-intronlen $min_intl --max-intronlen $max_intl --dta --phred33 -S $filename3.sam --met-file metrics.txt 2>> mapped.txt
            echo "samtools view -Sbo $filename3.bam $filename3.sam"
            samtools view -Sbo $filename3.bam $filename3.sam
            rm $filename3.sam
            mv mapped.txt metrics.txt "$bam_out"
            cufflinks_non_SRA
            mkdir index && mv $fbname* index
            echo "Cuffmerge only works with more than one file!!!"
            echo "##############################"
            echo "Pipeline executed successfully"
            echo "##############################"
            exit    
        elif [ $numb -gt 1 ] && [ "$cuff_merge" == 0 ]; then
            echo "------------------------------------------" >> mapped.txt
            echo "### Mapping percentages of" $filename3 "###" >> mapped.txt
            echo "------------------------------------------" >> mapped.txt
            echo "#######################"
            echo "Running Hisat2 mapping"
            echo "#######################"
            echo "hisat2 -x $fbname --rna-strandness $lib_type -1 ${filename}.fq.gz -2 ${filename2}.fq.gz -p $num_threads -5 $five_trim -3 $three_trim --min-intronlen $min_intl --max-intronlen $max_intl --dta --phred33 -S $filename3.sam --met-file metrics.txt 2>> mapped.txt"
            hisat2 -x $fbname --rna-strandness $lib_type -1 ${filename}.fq.gz -2 ${filename2}.fq.gz -p $num_threads -5 $five_trim -3 $three_trim --min-intronlen $min_intl --max-intronlen $max_intl --dta --phred33 -S $filename3.sam --met-file metrics.txt 2>> mapped.txt
            echo "samtools view -Sbo $filename3.bam $filename3.sam"
            samtools view -Sbo $filename3.bam $filename3.sam
            rm $filename3.sam
            cufflinks_non_SRA
        elif [ $numb -gt 1 ] && [ "$cuff_merge" != 0 ]; then
            echo "------------------------------------------" >> mapped.txt
            echo "### Mapping percentages of" $filename3 "###" >> mapped.txt
            echo "------------------------------------------" >> mapped.txt
            echo "#######################"
            echo "Running Hisat2 mapping"
            echo "#######################"
            echo "hisat2 -x $fbname --rna-strandness $lib_type -1 ${filename}.fq.gz -2 ${filename2}.fq.gz -p $num_threads -5 $five_trim -3 $three_trim --min-intronlen $min_intl --max-intronlen $max_intl --dta --phred33 -S $filename3.sam --met-file metrics.txt 2>> mapped.txt"
            hisat2 -x $fbname --rna-strandness $lib_type -1 ${filename}.fq.gz -2 ${filename2}.fq.gz -p $num_threads -5 $five_trim -3 $three_trim --min-intronlen $min_intl --max-intronlen $max_intl --dta --phred33 -S $filename3.sam --met-file metrics.txt 2>> mapped.txt
            echo "samtools view -Sbo $filename3.bam $filename3.sam"
            samtools view -Sbo $filename3.bam $filename3.sam
            rm $filename3.sam
            cufflinks_non_SRA
        fi

      elif [[ "$extension" =~ "fastq.gz" ]]; then
        filename=$(basename "$f" ".fastq.gz")
        filename2=${filename/_R1/_R2}
	      filename3=$(echo $filename | sed 's/_R1//')
        if [ $numb -eq 1 ] && [ "$cuff_merge" == 0 ]; then
            echo "------------------------------------------" >> mapped.txt
            echo "### Mapping percentages of" $filename3 "###" >> mapped.txt
            echo "------------------------------------------" >> mapped.txt
            echo "#######################"
            echo "Running Hisat2 mapping"
            echo "#######################"
            echo "hisat2 -x $fbname --rna-strandness $lib_type -1 ${filename}.fastq.gz -2 ${filename2}.fastq.gz -p $num_threads -5 $five_trim -3 $three_trim --min-intronlen $min_intl --max-intronlen $max_intl --dta --phred33 -S $filename3.sam --met-file metrics.txt 2>> mapped.txt"
            hisat2 -x $fbname --rna-strandness $lib_type -1 ${filename}.fastq.gz -2 ${filename2}.fastq.gz -p $num_threads -5 $five_trim -3 $three_trim --min-intronlen $min_intl --max-intronlen $max_intl --dta --phred33 -S $filename3.sam --met-file metrics.txt 2>> mapped.txt
            echo "samtools view -Sbo $filename3.bam $filename3.sam"
            samtools view -Sbo $filename3.bam $filename3.sam
            rm $filename3.sam
            mv mapped.txt metrics.txt "$bam_out"
            cufflinks_non_SRA
            mkdir index && mv $fbname* index
            echo "##############################"
            echo "Pipeline executed successfully"
            echo "##############################"
            exit
        elif [ $numb -eq 1 ] && [ "$cuff_merge" != 0 ]; then
            echo "------------------------------------------" >> mapped.txt
            echo "### Mapping percentages of" $filename3 "###" >> mapped.txt
            echo "------------------------------------------" >> mapped.txt
            echo "#######################"
            echo "Running Hisat2 mapping"
            echo "#######################"
            echo "hisat2 -x $fbname --rna-strandness $lib_type -1 ${filename}.fastq.gz -2 ${filename2}.fastq.gz -p $num_threads -5 $five_trim -3 $three_trim --min-intronlen $min_intl --max-intronlen $max_intl --dta --phred33 -S $filename3.sam --met-file metrics.txt 2>> mapped.txt"
            hisat2 -x $fbname --rna-strandness $lib_type -1 ${filename}.fastq.gz -2 ${filename2}.fastq.gz -p $num_threads -5 $five_trim -3 $three_trim --min-intronlen $min_intl --max-intronlen $max_intl --dta --phred33 -S $filename3.sam --met-file metrics.txt 2>> mapped.txt
            echo "samtools view -Sbo $filename3.bam $filename3.sam"
            samtools view -Sbo $filename3.bam $filename3.sam
            rm $filename3.sam
            mv mapped.txt metrics.txt "$bam_out"
            cufflinks_non_SRA
            mkdir index && mv $fbname* index
            echo "Cuffmerge only works with more than one file!!!"
            echo "##############################"
            echo "Pipeline executed successfully"
            echo "##############################"
            exit    
        elif [ $numb -gt 1 ] && [ "$cuff_merge" == 0 ]; then
            echo "------------------------------------------" >> mapped.txt
            echo "### Mapping percentages of" $filename3 "###" >> mapped.txt
            echo "------------------------------------------" >> mapped.txt
            echo "#######################"
            echo "Running Hisat2 mapping"
            echo "#######################"
            echo "hisat2 -x $fbname --rna-strandness $lib_type -1 ${filename}.fastq.gz -2 ${filename2}.fastq.gz -p $num_threads -5 $five_trim -3 $three_trim --min-intronlen $min_intl --max-intronlen $max_intl --dta --phred33 -S $filename3.sam --met-file metrics.txt 2>> mapped.txt"
            hisat2 -x $fbname --rna-strandness $lib_type -1 ${filename}.fastq.gz -2 ${filename2}.fastq.gz -p $num_threads -5 $five_trim -3 $three_trim --min-intronlen $min_intl --max-intronlen $max_intl --dta --phred33 -S $filename3.sam --met-file metrics.txt 2>> mapped.txt
            echo "samtools view -Sbo $filename3.bam $filename3.sam"
            samtools view -Sbo $filename3.bam $filename3.sam
            rm $filename3.sam
            cufflinks_non_SRA
        elif [ $numb -gt 1 ] && [ "$cuff_merge" != 0 ]; then
            echo "------------------------------------------" >> mapped.txt
            echo "### Mapping percentages of" $filename3 "###" >> mapped.txt
            echo "------------------------------------------" >> mapped.txt
            echo "#######################"
            echo "Running Hisat2 mapping"
            echo "#######################"
            echo "hisat2 -x $fbname --rna-strandness $lib_type -1 ${filename}.fastq.gz -2 ${filename2}.fastq.gz -p $num_threads -5 $five_trim -3 $three_trim --min-intronlen $min_intl --max-intronlen $max_intl --dta --phred33 -S $filename3.sam --met-file metrics.txt 2>> mapped.txt"
            hisat2 -x $fbname --rna-strandness $lib_type -1 ${filename}.fastq.gz -2 ${filename2}.fastq.gz -p $num_threads -5 $five_trim -3 $three_trim --min-intronlen $min_intl --max-intronlen $max_intl --dta --phred33 -S $filename3.sam --met-file metrics.txt 2>> mapped.txt
            echo "samtools view -Sbo $filename3.bam $filename3.sam"
            samtools view -Sbo $filename3.bam $filename3.sam
            rm $filename3.sam
            cufflinks_non_SRA
        fi

      elif [[ "$extension" =~ "fq" ]]; then
        filename=$(basename "$f" ".fq")
        filename2=${filename/_R1/_R2}
	      filename3=$(echo $filename | sed 's/_R1//')
        if [ $numb -eq 1 ] && [ "$cuff_merge" == 0 ]; then
            echo "------------------------------------------" >> mapped.txt
            echo "### Mapping percentages of" $filename3 "###" >> mapped.txt
            echo "------------------------------------------" >> mapped.txt
            echo "#######################"
            echo "Running Hisat2 mapping"
            echo "#######################"
            echo "hisat2 -x $fbname --rna-strandness $lib_type -1 ${filename}.fq -2 ${filename2}.fq -p $num_threads -5 $five_trim -3 $three_trim --min-intronlen $min_intl --max-intronlen $max_intl --dta --phred33 -S $filename3.sam --met-file metrics.txt 2>> mapped.txt"
            hisat2 -x $fbname --rna-strandness $lib_type -1 ${filename}.fq -2 ${filename2}.fq -p $num_threads -5 $five_trim -3 $three_trim --min-intronlen $min_intl --max-intronlen $max_intl --dta --phred33 -S $filename3.sam --met-file metrics.txt 2>> mapped.txt
            echo "samtools view -Sbo $filename3.bam $filename3.sam"
            samtools view -Sbo $filename3.bam $filename3.sam
            rm $filename3.sam
            mv mapped.txt metrics.txt "$bam_out"
            cufflinks_non_SRA
            mkdir index && mv $fbname* index
            echo "##############################"
            echo "Pipeline executed successfully"
            echo "##############################"
            exit
        elif [ $numb -eq 1 ] && [ "$cuff_merge" != 0 ]; then
            echo "------------------------------------------" >> mapped.txt
            echo "### Mapping percentages of" $filename3 "###" >> mapped.txt
            echo "------------------------------------------" >> mapped.txt
            echo "#######################"
            echo "Running Hisat2 mapping"
            echo "#######################"
            echo "hisat2 -x $fbname --rna-strandness $lib_type -1 ${filename}.fq -2 ${filename2}.fq -p $num_threads -5 $five_trim -3 $three_trim --min-intronlen $min_intl --max-intronlen $max_intl --dta --phred33 -S $filename3.sam --met-file metrics.txt 2>> mapped.txt"
            hisat2 -x $fbname --rna-strandness $lib_type -1 ${filename}.fq -2 ${filename2}.fq -p $num_threads -5 $five_trim -3 $three_trim --min-intronlen $min_intl --max-intronlen $max_intl --dta --phred33 -S $filename3.sam --met-file metrics.txt 2>> mapped.txt
            echo "samtools view -Sbo $filename3.bam $filename3.sam"
            samtools view -Sbo $filename3.bam $filename3.sam
            rm $filename3.sam
            mv mapped.txt metrics.txt "$bam_out"
            cufflinks_non_SRA
            mkdir index && mv $fbname* index
            echo "Cuffmerge only works with more than one file!!!"
            echo "##############################"
            echo "Pipeline executed successfully"
            echo "##############################"
            exit    
        elif [ $numb -gt 1 ] && [ "$cuff_merge" == 0 ]; then
            echo "------------------------------------------" >> mapped.txt
            echo "### Mapping percentages of" $filename3 "###" >> mapped.txt
            echo "------------------------------------------" >> mapped.txt
            echo "#######################"
            echo "Running Hisat2 mapping"
            echo "#######################"
            echo "hisat2 -x $fbname --rna-strandness $lib_type -1 ${filename}.fq -2 ${filename2}.fq -p $num_threads -5 $five_trim -3 $three_trim --min-intronlen $min_intl --max-intronlen $max_intl --dta --phred33 -S $filename3.sam --met-file metrics.txt 2>> mapped.txt"
            hisat2 -x $fbname --rna-strandness $lib_type -1 ${filename}.fq -2 ${filename2}.fq -p $num_threads -5 $five_trim -3 $three_trim --min-intronlen $min_intl --max-intronlen $max_intl --dta --phred33 -S $filename3.sam --met-file metrics.txt 2>> mapped.txt
            echo "samtools view -Sbo $filename3.bam $filename3.sam"
            samtools view -Sbo $filename3.bam $filename3.sam
            rm $filename3.sam
            cufflinks_non_SRA
        elif [ $numb -gt 1 ] && [ "$cuff_merge" != 0 ]; then
            echo "------------------------------------------" >> mapped.txt
            echo "### Mapping percentages of" $filename3 "###" >> mapped.txt
            echo "------------------------------------------" >> mapped.txt
            echo "#######################"
            echo "Running Hisat2 mapping"
            echo "#######################"
            echo "hisat2 -x $fbname --rna-strandness $lib_type -1 ${filename}.fq -2 ${filename2}.fq -p $num_threads -5 $five_trim -3 $three_trim --min-intronlen $min_intl --max-intronlen $max_intl --dta --phred33 -S $filename3.sam --met-file metrics.txt 2>> mapped.txt"
            hisat2 -x $fbname --rna-strandness $lib_type -1 ${filename}.fq -2 ${filename2}.fq -p $num_threads -5 $five_trim -3 $three_trim --min-intronlen $min_intl --max-intronlen $max_intl --dta --phred33 -S $filename3.sam --met-file metrics.txt 2>> mapped.txt
            echo "samtools view -Sbo $filename3.bam $filename3.sam"
            samtools view -Sbo $filename3.bam $filename3.sam
            rm $filename3.sam
            cufflinks_non_SRA
        fi

      elif [[ "$extension" =~ "fastq" ]]; then
        filename=$(basename "$f" ".fastq")
        filename2=${filename/_R1/_R2}
	      filename3=$(echo $filename | sed 's/_R1//')
        if [ $numb -eq 1 ] && [ "$cuff_merge" == 0 ]; then
            echo "------------------------------------------" >> mapped.txt
            echo "### Mapping percentages of" $filename3 "###" >> mapped.txt
            echo "------------------------------------------" >> mapped.txt
            echo "#######################"
            echo "Running Hisat2 mapping"
            echo "#######################"
            echo "hisat2 -x $fbname --rna-strandness $lib_type -1 ${filename}.fastq -2 ${filename2}.fastq -p $num_threads -5 $five_trim -3 $three_trim --min-intronlen $min_intl --max-intronlen $max_intl --dta --phred33 -S $filename3.sam --met-file metrics.txt 2>> mapped.txt"
            hisat2 -x $fbname --rna-strandness $lib_type -1 ${filename}.fastq -2 ${filename2}.fastq -p $num_threads -5 $five_trim -3 $three_trim --min-intronlen $min_intl --max-intronlen $max_intl --dta --phred33 -S $filename3.sam --met-file metrics.txt 2>> mapped.txt
            echo "samtools view -Sbo $filename3.bam $filename3.sam"
            samtools view -Sbo $filename3.bam $filename3.sam
            rm $filename3.sam
            mv mapped.txt metrics.txt "$bam_out"
            cufflinks_non_SRA
            mkdir index && mv $fbname* index
            echo "##############################"
            echo "Pipeline executed successfully"
            echo "##############################"
            exit
        elif [ $numb -eq 1 ] && [ "$cuff_merge" != 0 ]; then
            echo "------------------------------------------" >> mapped.txt
            echo "### Mapping percentages of" $filename3 "###" >> mapped.txt
            echo "------------------------------------------" >> mapped.txt
            echo "#######################"
            echo "Running Hisat2 mapping"
            echo "#######################"
            echo "hisat2 -x $fbname --rna-strandness $lib_type -1 ${filename}.fastq -2 ${filename2}.fastq -p $num_threads -5 $five_trim -3 $three_trim --min-intronlen $min_intl --max-intronlen $max_intl --dta --phred33 -S $filename3.sam --met-file metrics.txt 2>> mapped.txt"
            hisat2 -x $fbname --rna-strandness $lib_type -1 ${filename}.fastq -2 ${filename2}.fastq -p $num_threads -5 $five_trim -3 $three_trim --min-intronlen $min_intl --max-intronlen $max_intl --dta --phred33 -S $filename3.sam --met-file metrics.txt 2>> mapped.txt
            echo "samtools view -Sbo $filename3.bam $filename3.sam"
            samtools view -Sbo $filename3.bam $filename3.sam
            rm $filename3.sam
            mv mapped.txt metrics.txt "$bam_out"
            cufflinks_non_SRA
            mkdir index && mv $fbname* index
            echo "Cuffmerge only works with more than one file!!!"
            echo "##############################"
            echo "Pipeline executed successfully"
            echo "##############################"
            exit    
        elif [ $numb -gt 1 ] && [ "$cuff_merge" == 0 ]; then
            echo "------------------------------------------" >> mapped.txt
            echo "### Mapping percentages of" $filename3 "###" >> mapped.txt
            echo "------------------------------------------" >> mapped.txt
            echo "#######################"
            echo "Running Hisat2 mapping"
            echo "#######################"
            echo "hisat2 -x $fbname --rna-strandness $lib_type -1 ${filename}.fastq -2 ${filename2}.fastq -p $num_threads -5 $five_trim -3 $three_trim --min-intronlen $min_intl --max-intronlen $max_intl --dta --phred33 -S $filename3.sam --met-file metrics.txt 2>> mapped.txt"
            hisat2 -x $fbname --rna-strandness $lib_type -1 ${filename}.fastq -2 ${filename2}.fastq -p $num_threads -5 $five_trim -3 $three_trim --min-intronlen $min_intl --max-intronlen $max_intl --dta --phred33 -S $filename3.sam --met-file metrics.txt 2>> mapped.txt
            echo "samtools view -Sbo $filename3.bam $filename3.sam"
            samtools view -Sbo $filename3.bam $filename3.sam
            rm $filename3.sam
            cufflinks_non_SRA
        elif [ $numb -gt 1 ] && [ "$cuff_merge" != 0 ]; then
            echo "------------------------------------------" >> mapped.txt
            echo "### Mapping percentages of" $filename3 "###" >> mapped.txt
            echo "------------------------------------------" >> mapped.txt
            echo "#######################"
            echo "Running Hisat2 mapping"
            echo "#######################"
            echo "hisat2 -x $fbname --rna-strandness $lib_type -1 ${filename}.fastq -2 ${filename2}.fastq -p $num_threads -5 $five_trim -3 $three_trim --min-intronlen $min_intl --max-intronlen $max_intl --dta --phred33 -S $filename3.sam --met-file metrics.txt 2>> mapped.txt"
            hisat2 -x $fbname --rna-strandness $lib_type -1 ${filename}.fastq -2 ${filename2}.fastq -p $num_threads -5 $five_trim -3 $three_trim --min-intronlen $min_intl --max-intronlen $max_intl --dta --phred33 -S $filename3.sam --met-file metrics.txt 2>> mapped.txt
            echo "samtools view -Sbo $filename3.bam $filename3.sam"
            samtools view -Sbo $filename3.bam $filename3.sam
            rm $filename3.sam
            cufflinks_non_SRA
        fi   

      elif [ "$extension" != "fastq" ] || [ "$extension" != "fq" ] || [ "$extension" != "fastq.gz" ] || [ "$extension" != "fq.gz" ]; then
        echo "The extension" "$extension" "is not supported. Only .fq, .fq.gz, .fastq, .fastq.gz are only supported" 1>&2        
        exit 64
      fi 
    done
    
    if [ $numb -gt 1 ] && [ "$cuff_merge" == 0 ]; then
      mkdir index && mv $fbname* index
      mv mapped.txt metrics.txt "$bam_out"
      echo "##############################"
      echo "Pipeline executed successfully"
      echo "##############################"
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
            echo "------------------------------------------" >> mapped.txt
            echo "### Mapping percentages of" $filename3 "###" >> mapped.txt
            echo "------------------------------------------" >> mapped.txt
            echo "#######################"
            echo "Running Hisat2 mapping"
            echo "#######################"
            echo "hisat2 -x $fbname --rna-strandness $lib_type -1 ${filename}.fq.gz -2 ${filename2}.fq.gz -p $num_threads -5 $five_trim -3 $three_trim --min-intronlen $min_intl --max-intronlen $max_intl --dta --phred64 -S $filename3.sam --met-file metrics.txt 2>> mapped.txt"
            hisat2 -x $fbname --rna-strandness $lib_type -1 ${filename}.fq.gz -2 ${filename2}.fq.gz -p $num_threads -5 $five_trim -3 $three_trim --min-intronlen $min_intl --max-intronlen $max_intl --dta --phred64 -S $filename3.sam --met-file metrics.txt 2>> mapped.txt
            echo "samtools view -Sbo $filename3.bam $filename3.sam"
            samtools view -Sbo $filename3.bam $filename3.sam
            rm $filename3.sam
            mv mapped.txt metrics.txt "$bam_out"
            cufflinks_non_SRA
            mkdir index && mv $fbname* index
            echo "##############################"
            echo "Pipeline executed successfully"
            echo "##############################"
            exit
        elif [ $numb -eq 1 ] && [ "$cuff_merge" != 0 ]; then
            echo "------------------------------------------" >> mapped.txt
            echo "### Mapping percentages of" $filename3 "###" >> mapped.txt
            echo "------------------------------------------" >> mapped.txt
            echo "#######################"
            echo "Running Hisat2 mapping"
            echo "#######################"
            echo "hisat2 -x $fbname --rna-strandness $lib_type -1 ${filename}.fq.gz -2 ${filename2}.fq.gz -p $num_threads -5 $five_trim -3 $three_trim --min-intronlen $min_intl --max-intronlen $max_intl --dta --phred64 -S $filename3.sam --met-file metrics.txt 2>> mapped.txt"
            hisat2 -x $fbname --rna-strandness $lib_type -1 ${filename}.fq.gz -2 ${filename2}.fq.gz -p $num_threads -5 $five_trim -3 $three_trim --min-intronlen $min_intl --max-intronlen $max_intl --dta --phred64 -S $filename3.sam --met-file metrics.txt 2>> mapped.txt
            echo "samtools view -Sbo $filename3.bam $filename3.sam"
            samtools view -Sbo $filename3.bam $filename3.sam
            rm $filename3.sam
            mv mapped.txt metrics.txt "$bam_out"
            cufflinks_non_SRA
            mkdir index && mv $fbname* index
            echo "Cuffmerge only works with more than one file!!!"
            echo "##############################"
            echo "Pipeline executed successfully"
            echo "##############################"
            exit    
        elif [ $numb -gt 1 ] && [ "$cuff_merge" == 0 ]; then
            echo "------------------------------------------" >> mapped.txt
            echo "### Mapping percentages of" $filename3 "###" >> mapped.txt
            echo "------------------------------------------" >> mapped.txt
            echo "#######################"
            echo "Running Hisat2 mapping"
            echo "#######################"
            echo "hisat2 -x $fbname --rna-strandness $lib_type -1 ${filename}.fq.gz -2 ${filename2}.fq.gz -p $num_threads -5 $five_trim -3 $three_trim --min-intronlen $min_intl --max-intronlen $max_intl --dta --phred64 -S $filename3.sam --met-file metrics.txt 2>> mapped.txt"
            hisat2 -x $fbname --rna-strandness $lib_type -1 ${filename}.fq.gz -2 ${filename2}.fq.gz -p $num_threads -5 $five_trim -3 $three_trim --min-intronlen $min_intl --max-intronlen $max_intl --dta --phred64 -S $filename3.sam --met-file metrics.txt 2>> mapped.txt
            echo "samtools view -Sbo $filename3.bam $filename3.sam"
            samtools view -Sbo $filename3.bam $filename3.sam
            rm $filename3.sam
            cufflinks_non_SRA
        elif [ $numb -gt 1 ] && [ "$cuff_merge" != 0 ]; then
            echo "------------------------------------------" >> mapped.txt
            echo "### Mapping percentages of" $filename3 "###" >> mapped.txt
            echo "------------------------------------------" >> mapped.txt
            echo "#######################"
            echo "Running Hisat2 mapping"
            echo "#######################"
            echo "hisat2 -x $fbname --rna-strandness $lib_type -1 ${filename}.fq.gz -2 ${filename2}.fq.gz -p $num_threads -5 $five_trim -3 $three_trim --min-intronlen $min_intl --max-intronlen $max_intl --dta --phred64 -S $filename3.sam --met-file metrics.txt 2>> mapped.txt"
            hisat2 -x $fbname --rna-strandness $lib_type -1 ${filename}.fq.gz -2 ${filename2}.fq.gz -p $num_threads -5 $five_trim -3 $three_trim --min-intronlen $min_intl --max-intronlen $max_intl --dta --phred64 -S $filename3.sam --met-file metrics.txt 2>> mapped.txt
            echo "samtools view -Sbo $filename3.bam $filename3.sam"
            samtools view -Sbo $filename3.bam $filename3.sam
            rm $filename3.sam
            cufflinks_non_SRA
        fi

      elif [[ "$extension" =~ "fastq.gz" ]]; then
        filename=$(basename "$f" ".fastq.gz")
        filename2=${filename/_R1/_R2}
	      filename3=$(echo $filename | sed 's/_R1//')
        if [ $numb -eq 1 ] && [ "$cuff_merge" == 0 ]; then
            echo "------------------------------------------" >> mapped.txt
            echo "### Mapping percentages of" $filename3 "###" >> mapped.txt
            echo "------------------------------------------" >> mapped.txt
            echo "#######################"
            echo "Running Hisat2 mapping"
            echo "#######################"
            echo "hisat2 -x $fbname --rna-strandness $lib_type -1 ${filename}.fastq.gz -2 ${filename2}.fastq.gz -p $num_threads -5 $five_trim -3 $three_trim --min-intronlen $min_intl --max-intronlen $max_intl --dta --phred64 -S $filename3.sam --met-file metrics.txt 2>> mapped.txt"
            hisat2 -x $fbname --rna-strandness $lib_type -1 ${filename}.fastq.gz -2 ${filename2}.fastq.gz -p $num_threads -5 $five_trim -3 $three_trim --min-intronlen $min_intl --max-intronlen $max_intl --dta --phred64 -S $filename3.sam --met-file metrics.txt 2>> mapped.txt
            echo "samtools view -Sbo $filename3.bam $filename3.sam"
            samtools view -Sbo $filename3.bam $filename3.sam
            rm $filename3.sam
            mv mapped.txt metrics.txt "$bam_out"
            cufflinks_non_SRA
            mkdir index && mv $fbname* index
            echo "##############################"
            echo "Pipeline executed successfully"
            echo "##############################"
            exit
        elif [ $numb -eq 1 ] && [ "$cuff_merge" != 0 ]; then
            echo "------------------------------------------" >> mapped.txt
            echo "### Mapping percentages of" $filename3 "###" >> mapped.txt
            echo "------------------------------------------" >> mapped.txt
            echo "#######################"
            echo "Running Hisat2 mapping"
            echo "#######################"
            echo "hisat2 -x $fbname --rna-strandness $lib_type -1 ${filename}.fastq.gz -2 ${filename2}.fastq.gz -p $num_threads -5 $five_trim -3 $three_trim --min-intronlen $min_intl --max-intronlen $max_intl --dta --phred64 -S $filename3.sam --met-file metrics.txt 2>> mapped.txt"
            hisat2 -x $fbname --rna-strandness $lib_type -1 ${filename}.fastq.gz -2 ${filename2}.fastq.gz -p $num_threads -5 $five_trim -3 $three_trim --min-intronlen $min_intl --max-intronlen $max_intl --dta --phred64 -S $filename3.sam --met-file metrics.txt 2>> mapped.txt
            echo "samtools view -Sbo $filename3.bam $filename3.sam"
            samtools view -Sbo $filename3.bam $filename3.sam
            rm $filename3.sam
            mv mapped.txt metrics.txt "$bam_out"
            cufflinks_non_SRA
            mkdir index && mv $fbname* index
            echo "Cuffmerge only works with more than one file!!!"
            echo "##############################"
            echo "Pipeline executed successfully"
            echo "##############################"
            exit    
        elif [ $numb -gt 1 ] && [ "$cuff_merge" == 0 ]; then
            echo "------------------------------------------" >> mapped.txt
            echo "### Mapping percentages of" $filename3 "###" >> mapped.txt
            echo "------------------------------------------" >> mapped.txt
            echo "#######################"
            echo "Running Hisat2 mapping"
            echo "#######################"
            echo "hisat2 -x $fbname --rna-strandness $lib_type -1 ${filename}.fastq.gz -2 ${filename2}.fastq.gz -p $num_threads -5 $five_trim -3 $three_trim --min-intronlen $min_intl --max-intronlen $max_intl --dta --phred64 -S $filename3.sam --met-file metrics.txt 2>> mapped.txt"
            hisat2 -x $fbname --rna-strandness $lib_type -1 ${filename}.fastq.gz -2 ${filename2}.fastq.gz -p $num_threads -5 $five_trim -3 $three_trim --min-intronlen $min_intl --max-intronlen $max_intl --dta --phred64 -S $filename3.sam --met-file metrics.txt 2>> mapped.txt
            echo "samtools view -Sbo $filename3.bam $filename3.sam"
            samtools view -Sbo $filename3.bam $filename3.sam
            rm $filename3.sam
            cufflinks_non_SRA
        elif [ $numb -gt 1 ] && [ "$cuff_merge" != 0 ]; then
            echo "------------------------------------------" >> mapped.txt
            echo "### Mapping percentages of" $filename3 "###" >> mapped.txt
            echo "------------------------------------------" >> mapped.txt
            echo "#######################"
            echo "Running Hisat2 mapping"
            echo "#######################"
            echo "hisat2 -x $fbname --rna-strandness $lib_type -1 ${filename}.fastq.gz -2 ${filename2}.fastq.gz -p $num_threads -5 $five_trim -3 $three_trim --min-intronlen $min_intl --max-intronlen $max_intl --dta --phred64 -S $filename3.sam --met-file metrics.txt 2>> mapped.txt"
            hisat2 -x $fbname --rna-strandness $lib_type -1 ${filename}.fastq.gz -2 ${filename2}.fastq.gz -p $num_threads -5 $five_trim -3 $three_trim --min-intronlen $min_intl --max-intronlen $max_intl --dta --phred64 -S $filename3.sam --met-file metrics.txt 2>> mapped.txt
            echo "samtools view -Sbo $filename3.bam $filename3.sam"
            samtools view -Sbo $filename3.bam $filename3.sam
            rm $filename3.sam
            cufflinks_non_SRA
        fi

      elif [[ "$extension" =~ "fq" ]]; then
        filename=$(basename "$f" ".fq")
        filename2=${filename/_R1/_R2}
	      filename3=$(echo $filename | sed 's/_R1//')
        if [ $numb -eq 1 ] && [ "$cuff_merge" == 0 ]; then
            echo "------------------------------------------" >> mapped.txt
            echo "### Mapping percentages of" $filename3 "###" >> mapped.txt
            echo "------------------------------------------" >> mapped.txt
            echo "#######################"
            echo "Running Hisat2 mapping"
            echo "#######################"
            echo "hisat2 -x $fbname --rna-strandness $lib_type -1 ${filename}.fq -2 ${filename2}.fq -p $num_threads -5 $five_trim -3 $three_trim --min-intronlen $min_intl --max-intronlen $max_intl --dta --phred64 -S $filename3.sam --met-file metrics.txt 2>> mapped.txt"
            hisat2 -x $fbname --rna-strandness $lib_type -1 ${filename}.fq -2 ${filename2}.fq -p $num_threads -5 $five_trim -3 $three_trim --min-intronlen $min_intl --max-intronlen $max_intl --dta --phred64 -S $filename3.sam --met-file metrics.txt 2>> mapped.txt
            echo "samtools view -Sbo $filename3.bam $filename3.sam"
            samtools view -Sbo $filename3.bam $filename3.sam
            rm $filename3.sam
            mv mapped.txt metrics.txt "$bam_out"
            cufflinks_non_SRA
            mkdir index && mv $fbname* index
            echo "##############################"
            echo "Pipeline executed successfully"
            echo "##############################"
            exit
        elif [ $numb -eq 1 ] && [ "$cuff_merge" != 0 ]; then
            echo "------------------------------------------" >> mapped.txt
            echo "### Mapping percentages of" $filename3 "###" >> mapped.txt
            echo "------------------------------------------" >> mapped.txt
            echo "#######################"
            echo "Running Hisat2 mapping"
            echo "#######################"
            echo "hisat2 -x $fbname --rna-strandness $lib_type -1 ${filename}.fq -2 ${filename2}.fq -p $num_threads -5 $five_trim -3 $three_trim --min-intronlen $min_intl --max-intronlen $max_intl --dta --phred64 -S $filename3.sam --met-file metrics.txt 2>> mapped.txt"
            hisat2 -x $fbname --rna-strandness $lib_type -1 ${filename}.fq -2 ${filename2}.fq -p $num_threads -5 $five_trim -3 $three_trim --min-intronlen $min_intl --max-intronlen $max_intl --dta --phred64 -S $filename3.sam --met-file metrics.txt 2>> mapped.txt
            echo "samtools view -Sbo $filename3.bam $filename3.sam"
            samtools view -Sbo $filename3.bam $filename3.sam
            rm $filename3.sam
            mv mapped.txt metrics.txt "$bam_out"
            echo "Cuffmerge only works with more than one file!!!"
            cufflinks_non_SRA
            mkdir index && mv $fbname* index
            echo "##############################"
            echo "Pipeline executed successfully"
            echo "##############################"
            exit    
        elif [ $numb -gt 1 ] && [ "$cuff_merge" == 0 ]; then
            echo "------------------------------------------" >> mapped.txt
            echo "### Mapping percentages of" $filename3 "###" >> mapped.txt
            echo "------------------------------------------" >> mapped.txt
            echo "#######################"
            echo "Running Hisat2 mapping"
            echo "#######################"
            echo "hisat2 -x $fbname --rna-strandness $lib_type -1 ${filename}.fq -2 ${filename2}.fq -p $num_threads -5 $five_trim -3 $three_trim --min-intronlen $min_intl --max-intronlen $max_intl --dta --phred64 -S $filename3.sam --met-file metrics.txt 2>> mapped.txt"
            hisat2 -x $fbname --rna-strandness $lib_type -1 ${filename}.fq -2 ${filename2}.fq -p $num_threads -5 $five_trim -3 $three_trim --min-intronlen $min_intl --max-intronlen $max_intl --dta --phred64 -S $filename3.sam --met-file metrics.txt 2>> mapped.txt
            echo "samtools view -Sbo $filename3.bam $filename3.sam"
            samtools view -Sbo $filename3.bam $filename3.sam
            rm $filename3.sam
            cufflinks_non_SRA
        elif [ $numb -gt 1 ] && [ "$cuff_merge" != 0 ]; then
            echo "------------------------------------------" >> mapped.txt
            echo "### Mapping percentages of" $filename3 "###" >> mapped.txt
            echo "------------------------------------------" >> mapped.txt
            echo "#######################"
            echo "Running Hisat2 mapping"
            echo "#######################"
            echo "hisat2 -x $fbname --rna-strandness $lib_type -1 ${filename}.fq -2 ${filename2}.fq -p $num_threads -5 $five_trim -3 $three_trim --min-intronlen $min_intl --max-intronlen $max_intl --dta --phred64 -S $filename3.sam --met-file metrics.txt 2>> mapped.txt"
            hisat2 -x $fbname --rna-strandness $lib_type -1 ${filename}.fq -2 ${filename2}.fq -p $num_threads -5 $five_trim -3 $three_trim --min-intronlen $min_intl --max-intronlen $max_intl --dta --phred64 -S $filename3.sam --met-file metrics.txt 2>> mapped.txt
            echo "samtools view -Sbo $filename3.bam $filename3.sam"
            samtools view -Sbo $filename3.bam $filename3.sam
            rm $filename3.sam
            cufflinks_non_SRA
        fi

      elif [[ "$extension" =~ "fastq" ]]; then
        filename=$(basename "$f" ".fastq")
        filename2=${filename/_R1/_R2}
	      filename3=$(echo $filename | sed 's/_R1//')
        if [ $numb -eq 1 ] && [ "$cuff_merge" == 0 ]; then 
            echo "------------------------------------------" >> mapped.txt
            echo "### Mapping percentages of" $filename3 "###" >> mapped.txt
            echo "------------------------------------------" >> mapped.txt
            echo "#######################"
            echo "Running Hisat2 mapping"
            echo "#######################"
            echo "hisat2 -x $fbname --rna-strandness $lib_type -1 ${filename}.fastq -2 ${filename2}.fastq -p $num_threads -5 $five_trim -3 $three_trim --min-intronlen $min_intl --max-intronlen $max_intl --dta --phred64 -S $filename3.sam --met-file metrics.txt 2>> mapped.txt"
            hisat2 -x $fbname --rna-strandness $lib_type -1 ${filename}.fastq -2 ${filename2}.fastq -p $num_threads -5 $five_trim -3 $three_trim --min-intronlen $min_intl --max-intronlen $max_intl --dta --phred64 -S $filename3.sam --met-file metrics.txt 2>> mapped.txt
            echo "samtools view -Sbo $filename3.bam $filename3.sam"
            samtools view -Sbo $filename3.bam $filename3.sam
            rm $filename3.sam
            mv mapped.txt metrics.txt "$bam_out"
            cufflinks_non_SRA
            mkdir index && mv $fbname* index
            echo "##############################"
            echo "Pipeline executed successfully"
            echo "##############################"
            exit
        elif [ $numb -eq 1 ] && [ "$cuff_merge" != 0 ]; then
            echo "------------------------------------------" >> mapped.txt
            echo "### Mapping percentages of" $filename3 "###" >> mapped.txt
            echo "------------------------------------------" >> mapped.txt
            echo "#######################"
            echo "Running Hisat2 mapping"
            echo "#######################"
            echo "hisat2 -x $fbname --rna-strandness $lib_type -1 ${filename}.fastq -2 ${filename2}.fastq -p $num_threads -5 $five_trim -3 $three_trim --min-intronlen $min_intl --max-intronlen $max_intl --dta --phred64 -S $filename3.sam --met-file metrics.txt 2>> mapped.txt"
            hisat2 -x $fbname --rna-strandness $lib_type -1 ${filename}.fastq -2 ${filename2}.fastq -p $num_threads -5 $five_trim -3 $three_trim --min-intronlen $min_intl --max-intronlen $max_intl --dta --phred64 -S $filename3.sam --met-file metrics.txt 2>> mapped.txt
            echo "samtools view -Sbo $filename3.bam $filename3.sam"
            samtools view -Sbo $filename3.bam $filename3.sam
            rm $filename3.sam
            mv mapped.txt metrics.txt "$bam_out"
            mkdir index && mv $fbname* index
            cufflinks_non_SRA
            echo "Cuffmerge only works with more than one file!!!"
            echo "##############################"
            echo "Pipeline executed successfully"
            echo "##############################"
            exit    
        elif [ $numb -gt 1 ] && [ "$cuff_merge" == 0 ]; then
            echo "------------------------------------------" >> mapped.txt
            echo "### Mapping percentages of" $filename3 "###" >> mapped.txt
            echo "------------------------------------------" >> mapped.txt
            echo "#######################"
            echo "Running Hisat2 mapping"
            echo "#######################"
            echo "hisat2 -x $fbname --rna-strandness $lib_type -1 ${filename}.fastq -2 ${filename2}.fastq -p $num_threads -5 $five_trim -3 $three_trim --min-intronlen $min_intl --max-intronlen $max_intl --dta --phred64 -S $filename3.sam --met-file metrics.txt 2>> mapped.txt"
            hisat2 -x $fbname --rna-strandness $lib_type -1 ${filename}.fastq -2 ${filename2}.fastq -p $num_threads -5 $five_trim -3 $three_trim --min-intronlen $min_intl --max-intronlen $max_intl --dta --phred64 -S $filename3.sam --met-file metrics.txt 2>> mapped.txt
            echo "samtools view -Sbo $filename3.bam $filename3.sam"
            samtools view -Sbo $filename3.bam $filename3.sam
            rm $filename3.sam
            cufflinks_non_SRA
        elif [ $numb -gt 1 ] && [ "$cuff_merge" != 0 ]; then
            echo "------------------------------------------" >> mapped.txt
            echo "### Mapping percentages of" $filename3 "###" >> mapped.txt
            echo "------------------------------------------" >> mapped.txt
            echo "#######################"
            echo "Running Hisat2 mapping"
            echo "#######################"
            echo "hisat2 -x $fbname --rna-strandness $lib_type -1 ${filename}.fastq -2 ${filename2}.fastq -p $num_threads -5 $five_trim -3 $three_trim --min-intronlen $min_intl --max-intronlen $max_intl --dta --phred64 -S $filename3.sam --met-file metrics.txt 2>> mapped.txt"
            hisat2 -x $fbname --rna-strandness $lib_type -1 ${filename}.fastq -2 ${filename2}.fastq -p $num_threads -5 $five_trim -3 $three_trim --min-intronlen $min_intl --max-intronlen $max_intl --dta --phred64 -S $filename3.sam --met-file metrics.txt 2>> mapped.txt
            echo "samtools view -Sbo $filename3.bam $filename3.sam"
            samtools view -Sbo $filename3.bam $filename3.sam
            rm $filename3.sam
            cufflinks_non_SRA
        fi   

      elif [ "$extension" != "fastq" ] || [ "$extension" != "fq" ] || [ "$extension" != "fastq.gz" ] || [ "$extension" != "fq.gz" ]; then
        echo "The extension" "$extension" "is not supported. Only .fq, .fq.gz, .fastq, .fastq.gz are only supported" 1>&2        
        exit 64
      fi 
    done

    if [ $numb -gt 1 ] && [ "$cuff_merge" == 0 ]; then
      mkdir index && mv $fbname* index
      mv mapped.txt metrics.txt "$bam_out"
      echo "##############################"
      echo "Pipeline executed successfully"
      echo "##############################"
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
        echo "------------------------------------------" >> mapped.txt
        echo "### Mapping percentages of" $filename "###" >> mapped.txt
        echo "------------------------------------------" >> mapped.txt
        echo "#######################"
        echo "Running Hisat2 mapping"
        echo "#######################"
        echo "hisat2 -x $fbname --rna-strandness $lib_type -U $f -p $num_threads -5 $five_trim -3 $three_trim --min-intronlen $min_intl --max-intronlen $max_intl --dta --phred33 -S $filename.sam --met-file metrics.txt 2>> mapped.txt"
        hisat2 -x $fbname --rna-strandness $lib_type -U $f -p $num_threads -5 $five_trim -3 $three_trim --min-intronlen $min_intl --max-intronlen $max_intl --dta --phred33 -S $filename.sam --met-file metrics.txt 2>> mapped.txt
        echo "samtools view -Sbo $filename.bam $filename.sam"
        samtools view -Sbo $filename.bam $filename.sam
        rm $filename.sam
        mv mapped.txt metrics.txt "$bam_out"
        mkdir index && mv $fbname* index
        stringtie_non_SRA_single
        echo "##############################"
        echo "Pipeline executed successfully"
        echo "##############################"
        exit
      elif [ $numb -eq 1 ] && [ "$cuff_merge" != 0 ]; then
        extension=$(echo "$f" | sed -r 's/.*(fq|fq.gz|fastq|fastq.gz)$/\1/')
        filename=$(basename "$f" ".$extension")
        echo "------------------------------------------" >> mapped.txt
        echo "### Mapping percentages of" $filename "###" >> mapped.txt
        echo "------------------------------------------" >> mapped.txt
        echo "#######################"
        echo "Running Hisat2 mapping"
        echo "#######################"
        echo "hisat2 -x $fbname --rna-strandness $lib_type -U $f -p $num_threads -5 $five_trim -3 $three_trim --min-intronlen $min_intl --max-intronlen $max_intl --dta --phred33 -S $filename.sam --met-file metrics.txt 2>> mapped.txt"
        hisat2 -x $fbname --rna-strandness $lib_type -U $f -p $num_threads -5 $five_trim -3 $three_trim --min-intronlen $min_intl --max-intronlen $max_intl --dta --phred33 -S $filename.sam --met-file metrics.txt 2>> mapped.txt
        echo "samtools view -Sbo $filename.bam $filename.sam"
        samtools view -Sbo $filename.bam $filename.sam
        rm $filename.sam
        mv mapped.txt metrics.txt "$bam_out"
        mkdir index && mv $fbname* index
        stringtie_non_SRA_single
        echo "Cuffmerge only works with more than one file!!!"
        echo "##############################"
        echo "Pipeline executed successfully"
        echo "##############################"
        exit
      elif [ $numb -gt 1 ] && [ "$cuff_merge" == 0 ]; then
        extension=$(echo "$f" | sed -r 's/.*(fq|fq.gz|fastq|fastq.gz)$/\1/')
        filename=$(basename "$f" ".$extension")
        echo "------------------------------------------" >> mapped.txt
        echo "### Mapping percentages of" $filename "###" >> mapped.txt
        echo "------------------------------------------" >> mapped.txt
        echo "#######################"
        echo "Running Hisat2 mapping"
        echo "#######################"
        echo "hisat2 -x $fbname --rna-strandness $lib_type -U $f -p $num_threads -5 $five_trim -3 $three_trim --min-intronlen $min_intl --max-intronlen $max_intl --dta --phred33 -S $filename.sam --met-file metrics.txt 2>> mapped.txt"
        hisat2 -x $fbname --rna-strandness $lib_type -U $f -p $num_threads -5 $five_trim -3 $three_trim --min-intronlen $min_intl --max-intronlen $max_intl --dta --phred33 -S $filename.sam --met-file metrics.txt 2>> mapped.txt
        echo "samtools view -Sbo $filename.bam $filename.sam"
        samtools view -Sbo $filename.bam $filename.sam
        rm $filename.sam
        stringtie_non_SRA_single   
      elif [ $numb -gt 1 ] && [ "$cuff_merge" != 0 ]; then
        extension=$(echo "$f" | sed -r 's/.*(fq|fq.gz|fastq|fastq.gz)$/\1/')
        filename=$(basename "$f" ".$extension")
        echo "------------------------------------------" >> mapped.txt
        echo "### Mapping percentages of" $filename "###" >> mapped.txt
        echo "------------------------------------------" >> mapped.txt
        echo "#######################"
        echo "Running Hisat2 mapping"
        echo "#######################"
        echo "hisat2 -x $fbname --rna-strandness $lib_type -U $f -p $num_threads -5 $five_trim -3 $three_trim --min-intronlen $min_intl --max-intronlen $max_intl --dta --phred33 -S $filename.sam --met-file metrics.txt 2>> mapped.txt"
        hisat2 -x $fbname --rna-strandness $lib_type -U $f -p $num_threads -5 $five_trim -3 $three_trim --min-intronlen $min_intl --max-intronlen $max_intl --dta --phred33 -S $filename.sam --met-file metrics.txt 2>> mapped.txt
        echo "samtools view -Sbo $filename.bam $filename.sam"
        samtools view -Sbo $filename.bam $filename.sam
        rm $filename.sam
        stringtie_non_SRA_single
      fi 
    done

    if [ $numb -gt 1 ] && [ "$cuff_merge" == 0 ]; then
      mkdir index && mv $fbname* index
      mv mapped.txt metrics.txt "$bam_out"
      echo "##############################"
      echo "Pipeline executed successfully"
      echo "##############################"
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
        echo "------------------------------------------" >> mapped.txt
        echo "### Mapping percentages of" $filename "###" >> mapped.txt
        echo "------------------------------------------" >> mapped.txt
        echo "#######################"
        echo "Running Hisat2 mapping"
        echo "#######################"
        echo "hisat2 -x $fbname --rna-strandness $lib_type -U $f -p $num_threads -5 $five_trim -3 $three_trim --min-intronlen $min_intl --max-intronlen $max_intl --dta --phred64 -S $filename.sam --met-file metrics.txt 2>> mapped.txt"
        hisat2 -x $fbname --rna-strandness $lib_type -U $f -p $num_threads -5 $five_trim -3 $three_trim --min-intronlen $min_intl --max-intronlen $max_intl --dta --phred64 -S $filename.sam --met-file metrics.txt 2>> mapped.txt
        echo "samtools view -Sbo $filename.bam $filename.sam"
        samtools view -Sbo $filename.bam $filename.sam
        rm $filename.sam
        mv mapped.txt metrics.txt "$bam_out"
        mkdir index && mv $fbname* index
        stringtie_non_SRA_single
        echo "##############################"
        echo "Pipeline executed successfully"
        echo "##############################"
        exit 
      elif [ $numb -eq 1 ] && [ "$cuff_merge" != 0 ]; then
        extension=$(echo "$f" | sed -r 's/.*(fq|fq.gz|fastq|fastq.gz)$/\1/')
        filename=$(basename "$f" ".$extension")
        echo "------------------------------------------" >> mapped.txt
        echo "### Mapping percentages of" $filename "###" >> mapped.txt
        echo "------------------------------------------" >> mapped.txt
        echo "#######################"
        echo "Running Hisat2 mapping"
        echo "#######################"
        echo "hisat2 -x $fbname --rna-strandness $lib_type -U $f -p $num_threads -5 $five_trim -3 $three_trim --min-intronlen $min_intl --max-intronlen $max_intl --dta --phred64 -S $filename.sam --met-file metrics.txt 2>> mapped.txt"
        hisat2 -x $fbname --rna-strandness $lib_type -U $f -p $num_threads -5 $five_trim -3 $three_trim --min-intronlen $min_intl --max-intronlen $max_intl --dta --phred64 -S $filename.sam --met-file metrics.txt 2>> mapped.txt
        echo "samtools view -Sbo $filename.bam $filename.sam"
        samtools view -Sbo $filename.bam $filename.sam
        rm $filename.sam
        mv mapped.txt metrics.txt "$bam_out"
        mkdir index && mv $fbname* index
        stringtie_non_SRA_single
        echo "Cuffmerge only works with more than one file!!!"
        echo "##############################"
        echo "Pipeline executed successfully"
        echo "##############################"
        exit  
      elif [ $numb -gt 1 ] && [ "$cuff_merge" == 0 ]; then
        extension=$(echo "$f" | sed -r 's/.*(fq|fq.gz|fastq|fastq.gz)$/\1/')
        filename=$(basename "$f" ".$extension")
        echo "------------------------------------------" >> mapped.txt
        echo "### Mapping percentages of" $filename "###" >> mapped.txt
        echo "------------------------------------------" >> mapped.txt
        echo "#######################"
        echo "Running Hisat2 mapping"
        echo "#######################"
        echo "hisat2 -x $fbname --rna-strandness $lib_type -U $f -p $num_threads -5 $five_trim -3 $three_trim --min-intronlen $min_intl --max-intronlen $max_intl --dta --phred64 -S $filename.sam --met-file metrics.txt 2>> mapped.txt"
        hisat2 -x $fbname --rna-strandness $lib_type -U $f -p $num_threads -5 $five_trim -3 $three_trim --min-intronlen $min_intl --max-intronlen $max_intl --dta --phred64 -S $filename.sam --met-file metrics.txt 2>> mapped.txt
        echo "samtools view -Sbo $filename.bam $filename.sam"
        samtools view -Sbo $filename.bam $filename.sam
        rm $filename.sam
        stringtie_non_SRA_single
      elif [ $numb -gt 1 ] && [ "$cuff_merge" != 0 ]; then
        extension=$(echo "$f" | sed -r 's/.*(fq|fq.gz|fastq|fastq.gz)$/\1/')
        filename=$(basename "$f" ".$extension")
        echo "------------------------------------------" >> mapped.txt
        echo "### Mapping percentages of" $filename "###" >> mapped.txt
        echo "------------------------------------------" >> mapped.txt
        echo "#######################"
        echo "Running Hisat2 mapping"
        echo "#######################"
        echo "hisat2 -x $fbname --rna-strandness $lib_type -U $f -p $num_threads -5 $five_trim -3 $three_trim --min-intronlen $min_intl --max-intronlen $max_intl --dta --phred64 -S $filename.sam --met-file metrics.txt 2>> mapped.txt"
        hisat2 -x $fbname --rna-strandness $lib_type -U $f -p $num_threads -5 $five_trim -3 $three_trim --min-intronlen $min_intl --max-intronlen $max_intl --dta --phred64 -S $filename.sam --met-file metrics.txt 2>> mapped.txt
        echo "samtools view -Sbo $filename.bam $filename.sam"
        samtools view -Sbo $filename.bam $filename.sam
        rm $filename.sam
        stringtie_non_SRA_single
      fi 
    done
    
    if [ $numb -gt 1 ] && [ "$cuff_merge" == 0 ]; then
      mkdir index && mv $fbname* index
      mv mapped.txt metrics.txt "$bam_out"
      echo "##############################"
      echo "Pipeline executed successfully"
      echo "##############################"
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
        echo "------------------------------------------" >> mapped.txt
        echo "### Mapping percentages of" $filename "###" >> mapped.txt
        echo "------------------------------------------" >> mapped.txt
        echo "#######################"
        echo "Running Hisat2 mapping"
        echo "#######################"
        echo "hisat2 -x $fbname --rna-strandness $lib_type -U $f -p $num_threads -5 $five_trim -3 $three_trim --min-intronlen $min_intl --max-intronlen $max_intl --dta-cufflinks --phred33 -S $filename.sam --met-file metrics.txt 2>> mapped.txt"
        hisat2 -x $fbname --rna-strandness $lib_type -U $f -p $num_threads -5 $five_trim -3 $three_trim --min-intronlen $min_intl --max-intronlen $max_intl --dta-cufflinks --phred33 -S $filename.sam --met-file metrics.txt 2>> mapped.txt
        echo "samtools view -Sbo $filename.bam $filename.sam"
        samtools view -Sbo $filename.bam $filename.sam
        rm $filename.sam
        mv mapped.txt metrics.txt "$bam_out"
        mkdir index && mv $fbname* index
        cufflinks_non_SRA_single
        echo "##############################"
        echo "Pipeline executed successfully"
        echo "##############################"
        exit
      elif [ $numb -eq 1 ] && [ "$cuff_merge" != 0 ]; then 
        extension=$(echo "$f" | sed -r 's/.*(fq|fq.gz|fastq|fastq.gz)$/\1/')
        filename=$(basename "$f" ".$extension")
        echo "------------------------------------------" >> mapped.txt
        echo "### Mapping percentages of" $filename "###" >> mapped.txt
        echo "------------------------------------------" >> mapped.txt
        echo "#######################"
        echo "Running Hisat2 mapping"
        echo "#######################"
        echo "hisat2 -x $fbname --rna-strandness $lib_type -U $f -p $num_threads -5 $five_trim -3 $three_trim --min-intronlen $min_intl --max-intronlen $max_intl --dta-cufflinks --phred33 -S $filename.sam --met-file metrics.txt 2>> mapped.txt"
        hisat2 -x $fbname --rna-strandness $lib_type -U $f -p $num_threads -5 $five_trim -3 $three_trim --min-intronlen $min_intl --max-intronlen $max_intl --dta-cufflinks --phred33 -S $filename.sam --met-file metrics.txt 2>> mapped.txt
        echo "samtools view -Sbo $filename.bam $filename.sam"
        samtools view -Sbo $filename.bam $filename.sam
        rm $filename.sam
        mv mapped.txt metrics.txt "$bam_out"
        mkdir index && mv $fbname* index
        cufflinks_non_SRA_single
        echo "Cuffmerge only works with more than one file!!!"
        echo "##############################"
        echo "Pipeline executed successfully"
        echo "##############################"
        exit  
      elif [ $numb -gt 1 ] && [ "$cuff_merge" == 0 ]; then
        extension=$(echo "$f" | sed -r 's/.*(fq|fq.gz|fastq|fastq.gz)$/\1/')
        filename=$(basename "$f" ".$extension")
        echo "------------------------------------------" >> mapped.txt
        echo "### Mapping percentages of" $filename "###" >> mapped.txt
        echo "------------------------------------------" >> mapped.txt
        echo "#######################"
        echo "Running Hisat2 mapping"
        echo "#######################"
        echo "hisat2 -x $fbname --rna-strandness $lib_type -U $f -p $num_threads -5 $five_trim -3 $three_trim --min-intronlen $min_intl --max-intronlen $max_intl --dta-cufflinks --phred33 -S $filename.sam --met-file metrics.txt 2>> mapped.txt"
        hisat2 -x $fbname --rna-strandness $lib_type -U $f -p $num_threads -5 $five_trim -3 $three_trim --min-intronlen $min_intl --max-intronlen $max_intl --dta-cufflinks --phred33 -S $filename.sam --met-file metrics.txt 2>> mapped.txt
        echo "samtools view -Sbo $filename.bam $filename.sam"
        samtools view -Sbo $filename.bam $filename.sam
        rm $filename.sam
        cufflinks_non_SRA_single
      elif [ $numb -gt 1 ] && [ "$cuff_merge" != 0 ]; then
        extension=$(echo "$f" | sed -r 's/.*(fq|fq.gz|fastq|fastq.gz)$/\1/')
        filename=$(basename "$f" ".$extension")
        echo "------------------------------------------" >> mapped.txt
        echo "### Mapping percentages of" $filename "###" >> mapped.txt
        echo "------------------------------------------" >> mapped.txt
        echo "#######################"
        echo "Running Hisat2 mapping"
        echo "#######################"
        echo "hisat2 -x $fbname --rna-strandness $lib_type -U $f -p $num_threads -5 $five_trim -3 $three_trim --min-intronlen $min_intl --max-intronlen $max_intl --dta-cufflinks --phred33 -S $filename.sam --met-file metrics.txt 2>> mapped.txt"
        hisat2 -x $fbname --rna-strandness $lib_type -U $f -p $num_threads -5 $five_trim -3 $three_trim --min-intronlen $min_intl --max-intronlen $max_intl --dta-cufflinks --phred33 -S $filename.sam --met-file metrics.txt 2>> mapped.txt
        echo "samtools view -Sbo $filename.bam $filename.sam"
        samtools view -Sbo $filename.bam $filename.sam
        rm $filename.sam
        cufflinks_non_SRA_single
      fi          
    done
    
    if [ $numb -gt 1 ] && [ "$cuff_merge" == 0 ]; then
      mkdir index && mv $fbname* index
      mv mapped.txt metrics.txt "$bam_out"
      echo "##############################"
      echo "Pipeline executed successfully"
      echo "##############################"
      exit
    elif [ $numb -gt 1 ] && [ "$cuff_merge" != 0 ]; then
      cuff_merge_fun
    fi

# Phred 64

elif [ ! -z "$single_reads" ] && [ "$quality_64" != 0 ] && [ "$tra_as" == 0 ] && [ "$tra_cuff" != 0 ]; then
    mkdir "$bam_out"
    numb=$(ls "${single_reads[@]}" | wc -l)
    for f in "${single_reads[@]}"; do
      if [ $numb -eq 1 ] && [ "$cuff_merge" == 0 ]; then 
        extension=$(echo "$f" | sed -r 's/.*(fq|fq.gz|fastq|fastq.gz)$/\1/')
        filename=$(basename "$f" ".$extension")
        echo "------------------------------------------" >> mapped.txt
        echo "### Mapping percentages of" $filename "###" >> mapped.txt
        echo "------------------------------------------" >> mapped.txt
        echo "#######################"
        echo "Running Hisat2 mapping"
        echo "#######################"
        echo "hisat2 -x $fbname --rna-strandness $lib_type -U $f -p $num_threads -5 $five_trim -3 $three_trim --min-intronlen $min_intl --max-intronlen $max_intl --dta-cufflinks --phred64 -S $filename.sam --met-file metrics.txt 2>> mapped.txt"
        hisat2 -x $fbname --rna-strandness $lib_type -U $f -p $num_threads -5 $five_trim -3 $three_trim --min-intronlen $min_intl --max-intronlen $max_intl --dta-cufflinks --phred64 -S $filename.sam --met-file metrics.txt 2>> mapped.txt
        echo "samtools view -Sbo $filename.bam $filename.sam"
        samtools view -Sbo $filename.bam $filename.sam
        rm $filename.sam
        mv mapped.txt metrics.txt "$bam_out"
        mkdir index && mv $fbname* index
        cufflinks_non_SRA_single
        echo "##############################"
        echo "Pipeline executed successfully"
        echo "##############################"
        exit
      elif [ $numb -eq 1 ] && [ "$cuff_merge" != 0 ]; then 
        extension=$(echo "$f" | sed -r 's/.*(fq|fq.gz|fastq|fastq.gz)$/\1/')
        filename=$(basename "$f" ".$extension")
        echo "------------------------------------------" >> mapped.txt
        echo "### Mapping percentages of" $filename "###" >> mapped.txt
        echo "------------------------------------------" >> mapped.txt
        echo "#######################"
        echo "Running Hisat2 mapping"
        echo "#######################"
        echo "hisat2 -x $fbname --rna-strandness $lib_type -U $f -p $num_threads -5 $five_trim -3 $three_trim --min-intronlen $min_intl --max-intronlen $max_intl --dta-cufflinks --phred64 -S $filename.sam --met-file metrics.txt 2>> mapped.txt"
        hisat2 -x $fbname --rna-strandness $lib_type -U $f -p $num_threads -5 $five_trim -3 $three_trim --min-intronlen $min_intl --max-intronlen $max_intl --dta-cufflinks --phred64 -S $filename.sam --met-file metrics.txt 2>> mapped.txt
        echo "samtools view -Sbo $filename.bam $filename.sam"
        samtools view -Sbo $filename.bam $filename.sam
        rm $filename.sam
        mv mapped.txt metrics.txt "$bam_out"
        mkdir index && mv $fbname* index
        echo "Cuffmerge only works with more than one file!!!"
        cufflinks_non_SRA_single
        echo "##############################"
        echo "Pipeline executed successfully"
        echo "##############################"
        exit  
      elif [ $numb -gt 1 ] && [ "$cuff_merge" == 0 ]; then
        extension=$(echo "$f" | sed -r 's/.*(fq|fq.gz|fastq|fastq.gz)$/\1/')
        filename=$(basename "$f" ".$extension")
        echo "------------------------------------------" >> mapped.txt
        echo "### Mapping percentages of" $filename "###" >> mapped.txt
        echo "------------------------------------------" >> mapped.txt
        echo "#######################"
        echo "Running Hisat2 mapping"
        echo "#######################"
        echo "hisat2 -x $fbname --rna-strandness $lib_type -U $f -p $num_threads -5 $five_trim -3 $three_trim --min-intronlen $min_intl --max-intronlen $max_intl --dta-cufflinks --phred64 -S $filename.sam --met-file metrics.txt 2>> mapped.txt"
        hisat2 -x $fbname --rna-strandness $lib_type -U $f -p $num_threads -5 $five_trim -3 $three_trim --min-intronlen $min_intl --max-intronlen $max_intl --dta-cufflinks --phred64 -S $filename.sam --met-file metrics.txt 2>> mapped.txt
        echo "samtools view -Sbo $filename.bam $filename.sam"
        samtools view -Sbo $filename.bam $filename.sam
        rm $filename.sam
        cufflinks_non_SRA_single
      elif [ $numb -gt 1 ] && [ "$cuff_merge" != 0 ]; then
        extension=$(echo "$f" | sed -r 's/.*(fq|fq.gz|fastq|fastq.gz)$/\1/')
        filename=$(basename "$f" ".$extension")
        echo "------------------------------------------" >> mapped.txt
        echo "### Mapping percentages of" $filename "###" >> mapped.txt
        echo "------------------------------------------" >> mapped.txt
        echo "#######################"
        echo "Running Hisat2 mapping"
        echo "#######################"
        echo "hisat2 -x $fbname --rna-strandness $lib_type -U $f -p $num_threads -5 $five_trim -3 $three_trim --min-intronlen $min_intl --max-intronlen $max_intl --dta-cufflinks --phred64 -S $filename.sam --met-file metrics.txt 2>> mapped.txt"
        hisat2 -x $fbname --rna-strandness $lib_type -U $f -p $num_threads -5 $five_trim -3 $three_trim --min-intronlen $min_intl --max-intronlen $max_intl --dta-cufflinks --phred64 -S $filename.sam --met-file metrics.txt 2>> mapped.txt
        echo "samtools view -Sbo $filename.bam $filename.sam"
        samtools view -Sbo $filename.bam $filename.sam
        rm $filename.sam
        cufflinks_non_SRA_single
      fi          
    done

    if [ $numb -gt 1 ] && [ "$cuff_merge" == 0 ]; then
      mkdir index && mv $fbname* index
      mv mapped.txt metrics.txt "$bam_out"
      echo "##############################"
      echo "Pipeline executed successfully"
      echo "##############################"
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
      echo "------------------------------------------" >> mapped.txt
      echo "### Mapping percentages of" $f "###" >> mapped.txt
      echo "------------------------------------------" >> mapped.txt
      echo "#######################"
      echo "Running Hisat2 mapping"
      echo "#######################"
      echo "hisat2 -x $fbname --rna-strandness $lib_type --sra-acc $f -p $num_threads -5 $five_trim -3 $three_trim --min-intronlen $min_intl --max-intronlen $max_intl --dta --phred33 -S $f.sam --met-file metrics.txt 2>> mapped.txt"
      hisat2 -x $fbname --rna-strandness $lib_type --sra-acc $f -p $num_threads -5 $five_trim -3 $three_trim --min-intronlen $min_intl --max-intronlen $max_intl --dta --phred33 -S $f.sam --met-file metrics.txt 2>> mapped.txt
      echo "samtools view -Sbo $f.bam $f.sam"
      samtools view -Sbo $f.bam $f.sam
      rm $f.sam
      stringtie_SRA_multi      
      done < "$sra_id"
      if [ "$cuff_merge" == 0 ]; then
        mkdir index && mv $fbname* index
        mv mapped.txt metrics.txt "$bam_out"
        echo "##############################"
        echo "Pipeline executed successfully"
        echo "##############################"
        exit
      elif [ "$cuff_merge" != 0 ]; then
        cuff_merge_fun
      fi
  else    
    mkdir "$bam_out"
    echo "------------------------------------------" >> mapped.txt
    echo "### Mapping percentages of" $sra_id "###" >> mapped.txt
    echo "------------------------------------------" >> mapped.txt
    echo "#######################"
    echo "Running Hisat2 mapping"
    echo "#######################"
    echo "hisat2 -x $fbname --rna-strandness $lib_type --sra-acc $sra_id -p $num_threads -5 $five_trim -3 $three_trim --min-intronlen $min_intl --max-intronlen $max_intl --dta --phred33 -S $sra_id.sam --met-file metrics.txt 2>> mapped.txt"
    hisat2 -x $fbname --rna-strandness $lib_type --sra-acc $sra_id -p $num_threads -5 $five_trim -3 $three_trim --min-intronlen $min_intl --max-intronlen $max_intl --dta --phred33 -S $sra_id.sam --met-file metrics.txt 2>> mapped.txt
    echo "samtools view -Sbo $sra_id.bam $sra_id.sam"
    samtools view -Sbo $sra_id.bam $sra_id.sam
    rm $sra_id.sam
    mkdir index && mv $fbname* index
    mv mapped.txt metrics.txt "$bam_out"
    stringtie_SRA_single
    echo "##############################"
    echo "Pipeline executed successfully"
    echo "##############################"
    if [ "$cuff_merge" != 0 ]; then
      echo "cuffmerge only works with more than one SRA accesions. Use File containing SRA id's option to use cuffmerge"
    fi  
  fi

# phred 64

elif [ ! -z $sra_id ] && [ "$quality_64" != 0 ] && [ "$tra_as" != 0 ] && [ "$tra_cuff" == 0 ]; then
  if [[ -f $sra_id ]]; then
    mkdir "$bam_out"    
    echo "------------------------------------------" >> mapped.txt
    echo "### Mapping percentages of" $f "###" >> mapped.txt
    echo "------------------------------------------" >> mapped.txt
    echo "#######################"
    echo "Running Hisat2 mapping"
    echo "#######################"
    while read f; do
      echo "hisat2 -x $fbname --rna-strandness $lib_type --sra-acc $f -p $num_threads -5 $five_trim -3 $three_trim --min-intronlen $min_intl --max-intronlen $max_intl --dta --phred64 -S $f.sam --met-file metrics.txt 2>> mapped.txt"
      hisat2 -x $fbname --rna-strandness $lib_type --sra-acc $f -p $num_threads -5 $five_trim -3 $three_trim --min-intronlen $min_intl --max-intronlen $max_intl --dta --phred64 -S $f.sam --met-file metrics.txt 2>> mapped.txt
      echo "samtools view -Sbo $f.bam $f.sam"
      samtools view -Sbo $f.bam $f.sam
      rm $f.sam
      mv mapped.txt metrics.txt "$bam_out"
      stringtie_SRA_multi      
      done < "$sra_id"
      if [ "$cuff_merge" == 0 ]; then
        mkdir index && mv $fbname* index
        mv mapped.txt metrics.txt "$bam_out"
        echo "##############################"
        echo "Pipeline executed successfully"
        echo "##############################"
        exit
      elif [ "$cuff_merge" != 0 ]; then
        cuff_merge_fun
      fi
  else    
    mkdir "$bam_out"
    echo "------------------------------------------" >> mapped.txt
    echo "### Mapping percentages of" $sra_id "###" >> mapped.txt
    echo "------------------------------------------" >> mapped.txt
    echo "#######################"
    echo "Running Hisat2 mapping"
    echo "#######################"
    echo "hisat2 -x $fbname --rna-strandness $lib_type --sra-acc $sra_id -p $num_threads -5 $five_trim -3 $three_trim --min-intronlen $min_intl --max-intronlen $max_intl --dta --phred64 -S $sra_id.sam --met-file metrics.txt 2>> mapped.txt"
    hisat2 -x $fbname --rna-strandness $lib_type --sra-acc $sra_id -p $num_threads -5 $five_trim -3 $three_trim --min-intronlen $min_intl --max-intronlen $max_intl --dta --phred64 -S $sra_id.sam --met-file metrics.txt 2>> mapped.txt
    echo "samtools view -Sbo $sra_id.bam $sra_id.sam"
    samtools view -Sbo $sra_id.bam $sra_id.sam
    rm $sra_id.sam
    mkdir index && mv $fbname* index
    mv mapped.txt metrics.txt "$bam_out"
    stringtie_SRA_single
    echo "##############################"
    echo "Pipeline executed successfully"
    echo "##############################"
    if [ "$cuff_merge" != 0 ]; then
      echo "cuffmerge only works with more than one SRA accesions. Use File containing SRA id's option to use cuffmerge"
    fi
  fi 

# Cufflinks

# phred 33

elif [ ! -z $sra_id ] && [ "$quality_33" != 0 ] && [ "$tra_as" == 0 ] && [ "$tra_cuff" != 0 ]; then
  if [[ -f $sra_id ]]; then
    mkdir "$bam_out"
    echo "------------------------------------------" >> mapped.txt
    echo "### Mapping percentages of" $f "###" >> mapped.txt
    echo "------------------------------------------" >> mapped.txt
    echo "#######################"
    echo "Running Hisat2 mapping"
    echo "#######################"    
    while read f; do
      echo "hisat2 -x $fbname --rna-strandness $lib_type --sra-acc $f -p $num_threads -5 $five_trim -3 $three_trim --min-intronlen $min_intl --max-intronlen $max_intl --dta-cufflinks --phred33 -S $f.sam --met-file metrics.txt 2>> mapped.txt"
      hisat2 -x $fbname --rna-strandness $lib_type --sra-acc $f -p $num_threads -5 $five_trim -3 $three_trim --min-intronlen $min_intl --max-intronlen $max_intl --dta-cufflinks --phred33 -S $f.sam --met-file metrics.txt 2>> mapped.txt
      echo "samtools view -Sbo $f.bam $f.sam"
      samtools view -Sbo $f.bam $f.sam
      rm $f.sam
      mv mapped.txt metrics.txt "$bam_out"
      cufflinks_SRA_multi      
      done < "$sra_id"
      if [ "$cuff_merge" == 0 ]; then
        mkdir index && mv $fbname* index
        mv mapped.txt metrics.txt "$bam_out"
        echo "##############################"
        echo "Pipeline executed successfully"
        echo "##############################"
        exit
      elif [ "$cuff_merge" != 0 ]; then
        cuff_merge_fun
      fi
  else    
    mkdir "$bam_out"
    echo "------------------------------------------" >> mapped.txt
    echo "### Mapping percentages of" $sra_id "###" >> mapped.txt
    echo "------------------------------------------" >> mapped.txt
    echo "#######################"
    echo "Running Hisat2 mapping"
    echo "#######################"
    echo "hisat2 -x $fbname --rna-strandness $lib_type --sra-acc $sra_id -p $num_threads -5 $five_trim -3 $three_trim --min-intronlen $min_intl --max-intronlen $max_intl --dta-cufflinks --phred33 -S $sra_id.sam --met-file metrics.txt 2>> mapped.txt"
    hisat2 -x $fbname --rna-strandness $lib_type --sra-acc $sra_id -p $num_threads -5 $five_trim -3 $three_trim --min-intronlen $min_intl --max-intronlen $max_intl --dta-cufflinks --phred33 -S $sra_id.sam --met-file metrics.txt 2>> mapped.txt
    echo "samtools view -Sbo $sra_id.bam $sra_id.sam"
    samtools view -Sbo $sra_id.bam $sra_id.sam
    rm $sra_id.sam
    mkdir index && mv $fbname* index
    mv mapped.txt metrics.txt "$bam_out"
    cufflinks_SRA_single
    echo "##############################"
    echo "Pipeline executed successfully"
    echo "##############################"
    if [ "$cuff_merge" != 0 ]; then
      echo "cuffmerge only works with more than one SRA accesions. Use File containing SRA id's option to use cuffmerge"
    fi
  fi

# Phred 64

elif [ ! -z $sra_id ] && [ "$quality_64" != 0 ] && [ "$tra_as" == 0 ] && [ "$tra_cuff" != 0 ]; then
  if [[ -f $sra_id ]]; then
    mkdir "$bam_out"
    echo "------------------------------------------" >> mapped.txt
    echo "### Mapping percentages of" $f "###" >> mapped.txt
    echo "------------------------------------------" >> mapped.txt
    echo "#######################"
    echo "Running Hisat2 mapping"
    echo "#######################"    
    while read f; do
      echo "hisat2 -x $fbname --rna-strandness $lib_type --sra-acc $f -p $num_threads -5 $five_trim -3 $three_trim --min-intronlen $min_intl --max-intronlen $max_intl --dta-cufflinks --phred64 -S $f.sam --met-file metrics.txt 2>> mapped.txt"
      hisat2 -x $fbname --rna-strandness $lib_type --sra-acc $f -p $num_threads -5 $five_trim -3 $three_trim --min-intronlen $min_intl --max-intronlen $max_intl --dta-cufflinks --phred64 -S $f.sam --met-file metrics.txt 2>> mapped.txt
      echo "samtools view -Sbo $f.bam $f.sam"
      samtools view -Sbo $f.bam $f.sam
      rm $f.sam
      mv mapped.txt metrics.txt "$bam_out"
      cufflinks_SRA_multi      
      done < "$sra_id"
      if [ "$cuff_merge" == 0 ]; then
        mkdir index && mv $fbname* index
        mv mapped.txt metrics.txt "$bam_out"
        echo "##############################"
        echo "Pipeline executed successfully"
        echo "##############################"
        exit
      elif [ "$cuff_merge" != 0 ]; then
        cuff_merge_fun
      fi
  else    
    mkdir "$bam_out"
    echo "------------------------------------------" >> mapped.txt
    echo "### Mapping percentages of" $sra_id "###" >> mapped.txt
    echo "------------------------------------------" >> mapped.txt
    echo "#######################"
    echo "Running Hisat2 mapping"
    echo "#######################"
    echo "hisat2 -x $fbname --rna-strandness $lib_type --sra-acc $sra_id -p $num_threads -5 $five_trim -3 $three_trim --min-intronlen $min_intl --max-intronlen $max_intl --dta-cufflinks --phred64 -S $sra_id.sam --met-file metrics.txt 2>> mapped.txt"
    hisat2 -x $fbname --rna-strandness $lib_type --sra-acc $sra_id -p $num_threads -5 $five_trim -3 $three_trim --min-intronlen $min_intl --max-intronlen $max_intl --dta-cufflinks --phred64 -S $sra_id.sam --met-file metrics.txt 2>> mapped.txt
    echo "samtools view -Sbo $sra_id.bam $sra_id.sam"
    samtools view -Sbo $sra_id.bam $sra_id.sam
    rm $sra_id.sam
    mkdir index && mv $fbname* index
    mv mapped.txt metrics.txt "$bam_out"
    cufflinks_SRA_single
    echo "##############################"
    echo "Pipeline executed successfully"
    echo "##############################"
    if [ "$cuff_merge" != 0 ]; then
      echo "Cuffmerge only works with more than one SRA accesions. Use File containing SRA id's option to use cuffmerge!!!!"
    fi
  fi
fi

##### End ########
