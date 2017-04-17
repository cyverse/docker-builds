#!/bin/bash

usage() {
      echo ""
      echo "Usage : sh $0 {-g <reference_genome> | -u <custom_reference> | -i <Index_folder>} {-A <reference_annotation> | -R <custom reference annotation>} \
            -l lib_type {-1 <left_reads> -2 <right_reads> | -U <single_reads> | -s <sra_id>} -O <output_folder for Bam files> -p num_threads -5 <integer> \
            -3 <integer> [-q phred_33 -Q phred_64 -m min_intron -M max_intron -t transcriptome_stringtie -c cufflinks] -f <integer> -e <cuff_m>"
      echo ""

cat <<'EOF'

  -g <reference genome>

  -u </path/to/custom reference>

  -i </path/to/input folder>

  -A <refernce annotation>

  -R </path/to/cutsom reference annotation>

  -l Library type

  -1 </path/to/reads_1>

  -2 </path/to/reads_2>

  -U </path/to/single_reads>

  -O </path/to/bam output folder>

  -s SRA ID

  -p Number of threads
  
  -5 5' trim
 
  -3 3' trim

  -q phred33

  -Q phred64 

  -m Minimum intron length

  -M Maximum intron length

  -t StringTie (Report alignments tailored for transcript assemblers including StringTie)

  -c Cufflinks (Report alignments tailored specifically for Cufflinks)
  
  -f threshold

  -m Cuffmerge

EOF
    exit 0
}

quality_33=0
quality_64=0
tra_as=0
tra_cuff=0
cuff_merge=0


while getopts ":hg:u:i:A:R:l:1:2:U:O:s:p:5:3:f:qQtcem:M:" opt; do
  case $opt in
    g)
    referencegenome=$OPTARG
     ;;
    u)
    user_referencegenome=$OPTARG
     ;;  
    i)
    index_folder=$OPTARG # Input folder
     ;;
    A)
    referenceannotation=$OPTARG
     ;;
    R)
    user_referenceannotation=$OPTARG
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
    sra_id=$OPTARG # SRA ID
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
       rm listtoremove.txt
    elif [ "$threshold" -eq "$param4" ]; then
       grep " transcript" "$bam_out"/$sra_id.gtf | grep -e 'cov "3.' -e 'cov "2.' -e 'cov "1.' -e 'cov "0.' | cut -f 9 | cut -d " " -f 4 | sort -u > listtoremove.txt
       grep -vFf listtoremove.txt "$bam_out"/$sra_id.gtf >"$bam_out"/$sra_id.gtf.filtered.gtf
       rm listtoremove.txt
    elif [ "$threshold" -eq "$param3" ]; then
       grep " transcript" "$bam_out"/$sra_id.gtf | grep -e 'cov "2.' -e 'cov "1.' -e 'cov "0.' | cut -f 9 | cut -d " " -f 4 | sort -u > listtoremove.txt
       grep -vFf listtoremove.txt "$bam_out"/$sra_id.gtf >"$bam_out"/$sra_id.gtf.filtered.gtf
       rm listtoremove.txt
    elif [ "$threshold" -eq "$param2" ]; then
       grep " transcript" "$bam_out"/$sra_id.gtf | grep -e 'cov "1.' -e 'cov "0.' | cut -f 9 | cut -d " " -f 4 | sort -u > listtoremove.txt
       grep -vFf listtoremove.txt "$bam_out"/$sra_id.gtf >"$bam_out"/$sra_id.gtf.filtered.gtf
       rm listtoremove.txt
    elif [ "$threshold" -eq "$param1" ]; then
       grep " transcript" "$bam_out"/$sra_id.gtf | grep -e 'cov "0.' | cut -f 9 | cut -d " " -f 4 | sort -u > listtoremove.txt
       grep -vFf listtoremove.txt "$bam_out"/$sra_id.gtf >"$bam_out"/$sra_id.gtf.filtered.gtf
       rm listtoremove.txt
    elif [ "$threshold" -eq "$param0" ]; then
        mv "$bam_out"/$sra_id.gtf "$bam_out"/$sra_id.gtf.filtered.gtf
    else
        echo "Invalid coverage parameter. Please select a whole number between 0-5"
        exit
    fi
}

stringtie_SRA_single()
{
    samtools sort -O BAM -T temp_files $sra_id.sam -o "$bam_out"/$sra_id.sorted.bam --threads $num_threads
    rm $sra_id.sam
    if [ ! -z "$user_referenceannotation" ] && [ -z "$referenceannotation" ]; then
      stringtie -G $user_referenceannotation "$bam_out"/$sra_id.sorted.bam -o "$bam_out"/$sra_id.gtf -p $num_threads
      coverge_cuffoff_SRA_single  
      cuffcompare "$bam_out"/$sra_id.gtf.filtered.gtf -r $user_referenceannotation -o $sra_id
      mv *.tracking *.loci *.combined.gtf *.stats "$bam_out" 
    elif [ -z "$user_referenceannotation" ] && [ ! -z "$referenceannotation" ]; then
      stringtie -G $referenceannotation "$bam_out"/$sra_id.sorted.bam -o "$bam_out"/$sra_id.gtf -p $num_threads
      coverge_cuffoff_SRA_single  
      cuffcompare "$bam_out"/$sra_id.gtf.filtered.gtf -r $referenceannotation -o $sra_id
      mv *.tracking *.loci *.combined.gtf *.stats "$bam_out"
    fi
}

cufflinks_SRA_single()
{
    samtools sort -O BAM -T temp_files $sra_id.sam -o "$bam_out"/$sra_id.sorted.bam --threads $num_threads
    rm $sra_id.sam
    if [ ! -z "$user_referenceannotation" ] && [ -z "$referenceannotation" ]; then
      cufflinks "$bam_out"/$sra_id.sorted.bam -p $num_threads -g $user_referenceannotation -o "$bam_out"
      mv "$bam_out"/transcripts.gtf "$bam_out"/$sra_id.gtf
      mv "$bam_out"/isoforms.fpkm_tracking "$bam_out"/$sra_id.isoforms.fpkm_tracking
      mv "$bam_out"/genes.fpkm_tracking "$bam_out"/$sra_id.genes.fpkm_tracking
      coverge_cuffoff_SRA_single        
      cuffcompare "$bam_out"/$sra_id.gtf.filtered.gtf -r $user_referenceannotation -o $sra_id
      mv *.tracking *.loci *.combined.gtf *.stats "$bam_out" 
    elif [ -z "$user_referenceannotation" ] && [ ! -z "$referenceannotation" ]; then
      cufflinks "$bam_out"/$sra_id.sorted.bam -p $num_threads -g $referenceannotation -o "$bam_out"
      mv "$bam_out"/transcripts.gtf "$bam_out"/$sra_id.gtf
      mv "$bam_out"/isoforms.fpkm_tracking "$bam_out"/$sra_id.isoforms.fpkm_tracking
      mv "$bam_out"/genes.fpkm_tracking "$bam_out"/$sra_id.genes.fpkm_tracking
      coverge_cuffoff_SRA_single           
      cuffcompare "$bam_out"/$sra_id.gtf.filtered.gtf -r $user_referenceannotation -o $sra_id
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
       rm listtoremove.txt
    elif [ "$threshold" -eq "$param4" ]; then
       grep " transcript" "$bam_out"/$f.gtf | grep -e 'cov "3.' -e 'cov "2.' -e 'cov "1.' -e 'cov "0.' | cut -f 9 | cut -d " " -f 4 | sort -u > listtoremove.txt
       grep -vFf listtoremove.txt "$bam_out"/$f.gtf >"$bam_out"/$f.gtf.filtered.gtf
       rm listtoremove.txt
    elif [ "$threshold" -eq "$param3" ]; then
       grep " transcript" "$bam_out"/$f.gtf | grep -e 'cov "2.' -e 'cov "1.' -e 'cov "0.' | cut -f 9 | cut -d " " -f 4 | sort -u > listtoremove.txt
       grep -vFf listtoremove.txt "$bam_out"/$f.gtf >"$bam_out"/$f.gtf.filtered.gtf
       rm listtoremove.txt
    elif [ "$threshold" -eq "$param2" ]; then
       grep " transcript" "$bam_out"/$f.gtf | grep -e 'cov "1.' -e 'cov "0.' | cut -f 9 | cut -d " " -f 4 | sort -u > listtoremove.txt
       grep -vFf listtoremove.txt "$bam_out"/$f.gtf >"$bam_out"/$f.gtf.filtered.gtf
       rm listtoremove.txt
    elif [ "$threshold" -eq "$param1" ]; then
       grep " transcript" "$bam_out"/$f.gtf | grep -e 'cov "0.' | cut -f 9 | cut -d " " -f 4 | sort -u > listtoremove.txt
       grep -vFf listtoremove.txt "$bam_out"/$f.gtf >"$bam_out"/$f.gtf.filtered.gtf
       rm listtoremove.txt
    elif [ "$threshold" -eq "$param0" ]; then
        mv "$bam_out"/$f.gtf "$bam_out"/$f.gtf.filtered.gtf
    else
        echo "Invalid coverage parameter. Please select a whole number between 0-5"
        exit
    fi
}

stringtie_SRA_multi()
{
    samtools sort -O BAM -T temp_files $f.sam -o "$bam_out"/$f.sorted.bam --threads $num_threads
    rm $f.sam
    if [ ! -z "$user_referenceannotation" ] && [ -z "$referenceannotation" ]; then
      stringtie -G $user_referenceannotation "$bam_out"/$f.sorted.bam -o "$bam_out"/$f.gtf -p $num_threads
      coverge_cuffoff_SRA_multi  
      cuffcompare "$bam_out"/$f.gtf.filtered.gtf -r $user_referenceannotation -o $f
      mv *.tracking *.loci *.combined.gtf *.stats "$bam_out"
    elif [ -z "$user_referenceannotation" ] && [ ! -z "$referenceannotation" ]; then
      stringtie -G $referenceannotation "$bam_out"/$f.sorted.bam -o "$bam_out"/$f.gtf -p $num_threads
      coverge_cuffoff_SRA_multi  
      cuffcompare "$bam_out"/$f.gtf.filtered.gtf -r $referenceannotation -o $f
      mv *.tracking *.loci *.combined.gtf *.stats "$bam_out"
    fi
}

cufflinks_SRA_multi()
{
    samtools sort -O BAM -T temp_files $f.sam -o "$bam_out"/$f.sorted.bam --threads $num_threads
    rm $f.sam
    if [ ! -z "$user_referenceannotation" ] && [ -z "$referenceannotation" ]; then
        cufflinks "$bam_out"/$f.sorted.bam -p $num_threads -g $user_referenceannotation -o "$bam_out"
        mv "$bam_out"/transcripts.gtf "$bam_out"/$f.gtf
        mv "$bam_out"/isoforms.fpkm_tracking "$bam_out"/$f.isoforms.fpkm_tracking
        mv "$bam_out"/skipped.gtf  "$bam_out"/$f.skipped.gtf
        coverge_cuffoff_SRA_multi  
        cuffcompare "$bam_out"/$f.gtf.filtered.gtf -r $user_referenceannotation -o $f
        mv *.tracking *.loci *.combined.gtf "$bam_out"
    elif [ -z "$user_referenceannotation" ] && [ ! -z "$referenceannotation" ]; then
        cufflinks "$bam_out"/$f.sorted.bam -p $num_threads -g $referenceannotation -o "$bam_out"
        mv "$bam_out"/transcripts.gtf "$bam_out"/$f.gtf
        mv "$bam_out"/isoforms.fpkm_tracking "$bam_out"/$f.isoforms.fpkm_tracking
        mv "$bam_out"/skipped.gtf  "$bam_out"/$f.skipped.gtf
        coverge_cuffoff_SRA_multi  
        cuffcompare "$bam_out"/$f.gtf.filtered.gtf -r $user_referenceannotation -o $f
        mv *.tracking *.loci *.combined.gtf *.stats "$bam_out"  
    fi
}

##################
### No SRA #######
##################

coverge_cuffoff_non_SRA()
{
    if [ "$threshold" -eq "$param5" ]; then
       grep " transcript" "$bam_out"/$filename.gtf | grep -e 'cov "4.' -e 'cov "3.' -e 'cov "2.' -e 'cov "1.' -e 'cov "0.' | cut -f 9 | cut -d " " -f 4 | sort -u > listtoremove.txt
       grep -vFf listtoremove.txt "$bam_out"/$filename.gtf >"$bam_out"/$filename.gtf.filtered.gtf
       rm listtoremove.txt
    elif [ "$threshold" -eq "$param4" ]; then
       grep " transcript" "$bam_out"/$filename.gtf | grep -e 'cov "3.' -e 'cov "2.' -e 'cov "1.' -e 'cov "0.' | cut -f 9 | cut -d " " -f 4 | sort -u > listtoremove.txt
       grep -vFf listtoremove.txt "$bam_out"/$filename.gtf >"$bam_out"/$filename.gtf.filtered.gtf
       rm listtoremove.txt
    elif [ "$threshold" -eq "$param3" ]; then
       grep " transcript" "$bam_out"/$filename.gtf | grep -e 'cov "2.' -e 'cov "1.' -e 'cov "0.' | cut -f 9 | cut -d " " -f 4 | sort -u > listtoremove.txt
       grep -vFf listtoremove.txt "$bam_out"/$filename.gtf >"$bam_out"/$filename.gtf.filtered.gtf
       rm listtoremove.txt
    elif [ "$threshold" -eq "$param2" ]; then
       grep " transcript" "$bam_out"/$filename.gtf | grep -e 'cov "1.' -e 'cov "0.' | cut -f 9 | cut -d " " -f 4 | sort -u > listtoremove.txt
       grep -vFf listtoremove.txt "$bam_out"/$filename.gtf >"$bam_out"/$filename.gtf.filtered.gtf
       rm listtoremove.txt
    elif [ "$threshold" -eq "$param1" ]; then
       grep " transcript" "$bam_out"/$filename.gtf | grep -e 'cov "0.' | cut -f 9 | cut -d " " -f 4 | sort -u > listtoremove.txt
       grep -vFf listtoremove.txt "$bam_out"/$filename.gtf >"$bam_out"/$filename.gtf.filtered.gtf
       rm listtoremove.txt
    elif [ "$threshold" -eq "$param0" ]; then
        mv "$bam_out"/$filename.gtf "$bam_out"/$filename.gtf.filtered.gtf
    else
        echo "Invalid coverage parameter. Please select a whole number between 0-5"
        exit
    fi
}

stringtie_non_SRA() 
{
    samtools sort -O BAM -T temp_files $filename.sam -o "$bam_out"/$filename.sorted.bam --threads $num_threads
    rm $filename.sam
    if [ ! -z "$user_referenceannotation" ] && [ -z "$referenceannotation" ]; then
        stringtie -G $user_referenceannotation "$bam_out"/$filename.sorted.bam -o "$bam_out"/$filename.gtf -p $num_threads
        coverge_cuffoff_non_SRA
        cuffcompare "$bam_out"/$filename.gtf.filtered.gtf -r $user_referenceannotation -o $filename
        mv *.tracking *.loci *.combined.gtf *.stats "$bam_out"
    elif [ -z "$user_referenceannotation" ] && [ ! -z "$referenceannotation" ]; then
        stringtie -G $referenceannotation "$bam_out"/$filename.sorted.bam -o "$bam_out"/$filename.gtf -p $num_threads
        coverge_cuffoff_non_SRA 
        cuffcompare "$bam_out"/$filename.gtf.filtered.gtf -r $referenceannotation -o $filename
        mv *.tracking *.loci *.combined.gtf *.stats "$bam_out"
    fi
}

cufflinks_non_SRA()
{
    samtools sort -O BAM -T temp_files $filename.sam -o "$bam_out"/$filename.sorted.bam --threads $num_threads
    rm $filename.sam
    if [ ! -z "$user_referenceannotation" ] && [ -z "$referenceannotation" ]; then
        cufflinks "$bam_out"/$filename.sorted.bam -p $num_threads -g $user_referenceannotation -o "$bam_out"
        mv "$bam_out"/skipped.gtf "$bam_out"/$filename.skipped.gtf
        mv "$bam_out"/transcripts.gtf "$bam_out"/$filename.gtf
        mv "$bam_out"/isoforms.fpkm_tracking "$bam_out"/$filename.isoforms.fpkm_tracking
        mv "$bam_out"/genes.fpkm_tracking "$bam_out"/$filename.genes.fpkm_tracking
        coverge_cuffoff_non_SRA        
        cuffcompare "$bam_out"/$filename.gtf.filtered.gtf -r $user_referenceannotation -o $filename
        mv *.tracking *.loci *.combined.gtf "$bam_out"
    elif [ -z "$user_referenceannotation" ] && [ ! -z "$referenceannotation" ]; then
        cufflinks "$bam_out"/$filename.sorted.bam -p $num_threads -g $referenceannotation -o "$bam_out"
        mv "$bam_out"/skipped.gtf "$bam_out"/$filename.skipped.gtf
        mv "$bam_out"/transcripts.gtf "$bam_out"/$filename.gtf
        mv "$bam_out"/isoforms.fpkm_tracking "$bam_out"/$filename.isoforms.fpkm_tracking
        mv "$bam_out"/genes.fpkm_tracking "$bam_out"/$filename.genes.fpkm_tracking
        coverge_cuffoff_non_SRA        
        cuffcompare "$bam_out"/$filename.gtf.filtered.gtf -r $user_referenceannotation -o $filename
        mv *.tracking *.loci *.combined.gtf *.stats "$bam_out"
    fi  
}

###############
## cuffmerge ##
###############

cuff_merge_fun()
{
    if [ ! -z "$user_referenceannotation" ] && [ -z "$referenceannotation" ] && [ "$cuff_merge" != 0 ]; then
        ls "$bam_out"/*combined.gtf | tr "\t" "\n" >> "$bam_out"/gtf_file.txt 
        cuffmerge -o "$bam_out"/merged_out -g $user_referenceannotation "$bam_out"/gtf_file.txt -p $num_threads
    elif [ -z "$user_referenceannotation" ] && [ ! -z "$referenceannotation" ] && [ "$cuff_merge" != 0 ]; then
        ls "$bam_out"/*combined.gtf | tr "\t" "\n" >> "$bam_out"/gtf_file.txt 
        cuffmerge -o "$bam_out"/merged_out -g $referenceannotation "$bam_out"/gtf_file.txt -p $num_threads
    fi
}

# ############################################################################################################################################################################################################################
# # Reference genome/Index
# ############################################################################################################################################################################################################################
# Index folder

if [ ! -z "$referencegenome" ] && [ -z "$user_referencegenome" ] && [ -z "$index_folder" ]; then
  hisat2-build $referencegenome ref_genome -p $num_threads
  fbname=$(basename "ref_genome" .ht2 | cut -d. -f1)
  $fbname

# Custom reference
elif [ -z "$referencegenome" ] && [ ! -z "$user_referencegenome" ] && [ -z "$index_folder" ]; then
  cp "$user_referencegenome" .
  hisat2-build genome.fas ref_genome -p $num_threads
  fbname=$(basename "ref_genome" .ht2 | cut -d. -f1)
  rm genome.fas

# Reference genomes from CyVerse
elif [ -z "$referencegenome" ] && [ -z "$user_referencegenome" ] && [ ! -z "$index_folder" ]; then
  for i in $index_folder/*; do
      cp $i .
     fbname=$(basename "$i" .ht2 | cut -d. -f1)
  done
fi

# ############################################################################################################################################################################################################################
# # Paired end reads
# ############################################################################################################################################################################################################################

# Phred 33

# Stringtie

if [ ! -z "$left_reads" ] && [ ! -z "$right_reads" ] && [ "$quality_33" != 0 ] && [ "$tra_as" != 0 ] && [ "$tra_cuff" == 0 ]; then
    mkdir "$bam_out"
    numb=$(ls "${left_reads[@]}" | wc -l)
    for f in "${left_reads[@]}"; do
      extension=${f#*.}
      if [ "$extension" == "fastq" ]; then
        if [ $numb -eq 1 ] && [ "$cuff_merge" == 0 ]; then
            filename=$(basename "$f" "_R1_001.fastq")
            hisat2 -x $fbname --rna-strandness $lib_type -1 ${filename}_R1_001.fastq -2 ${filename}_R2_001.fastq -S $filename.sam -p $num_threads -5 $five_trim -3 $three_trim --min-intronlen $min_intl --max-intronlen $max_intl --dta --phred33
            stringtie_non_SRA
            rm $fbname*
            exit
        elif [ $numb -eq 1 ] && [ "$cuff_merge" != 0 ]; then
            filename=$(basename "$f" "_R1_001.fastq")
            hisat2 -x $fbname --rna-strandness $lib_type -1 ${filename}_R1_001.fastq -2 ${filename}_R2_001.fastq -S $filename.sam -p $num_threads -5 $five_trim -3 $three_trim --min-intronlen $min_intl --max-intronlen $max_intl --dta --phred33
            stringtie_non_SRA
            rm $fbname*
            echo "cuffmerge only works with more than one file"
            exit    
        elif [ $numb -gt 1 ] && [ "$cuff_merge" != 0 ]; then
            filename=$(basename "$f" "_R1_001.fastq")
            hisat2 -x $fbname --rna-strandness $lib_type -1 ${filename}_R1_001.fastq -2 ${filename}_R2_001.fastq -S $filename.sam -p $num_threads -5 $five_trim -3 $three_trim --min-intronlen $min_intl --max-intronlen $max_intl --dta --phred33
            stringtie_non_SRA
        fi        
      elif [ "$extension" == "fq" ]; then
        if [ $numb -eq 1 ] && [ "$cuff_merge" == 0 ]; then
            filename=$(basename "$f" "_R1_001.fq")
            hisat2 -x $fbname --rna-strandness $lib_type -1 ${filename}_R1_001.fq -2 ${filename}_R2_001.fq -S $filename.sam -p $num_threads -5 $five_trim -3 $three_trim --min-intronlen $min_intl --max-intronlen $max_intl --dta --phred33
            stringtie_non_SRA
            rm $fbname*
            exit
        elif [ $numb -eq 1 ] && [ "$cuff_merge" != 0 ]; then
            filename=$(basename "$f" "_R1_001.fq")
            hisat2 -x $fbname --rna-strandness $lib_type -1 ${filename}_R1_001.fq -2 ${filename}_R2_001.fq -S $filename.sam -p $num_threads -5 $five_trim -3 $three_trim --min-intronlen $min_intl --max-intronlen $max_intl --dta --phred33
            stringtie_non_SRA
            rm $fbname*
            echo "cuffmerge only works with more than one file"
            exit    
        elif [ $numb -gt 1 ] && [ "$cuff_merge" != 0 ]; then
            filename=$(basename "$f" "_R1_001.fq")
            hisat2 -x $fbname --rna-strandness $lib_type -1 ${filename}_R1_001.fq -2 ${filename}_R2_001.fq -S $filename.sam -p $num_threads -5 $five_trim -3 $three_trim --min-intronlen $min_intl --max-intronlen $max_intl --dta --phred33
            stringtie_non_SRA
        fi       
      elif [ "$extension" == "fastq.gz" ]; then
        if [ $numb -eq 1 ] && [ "$cuff_merge" == 0 ]; then
            filename=$(basename "$f" "_R1_001.fastq.gz")
            hisat2 -x $fbname --rna-strandness $lib_type -1 ${filename}_R1_001.fastq.gz -2 ${filename}_R2_001.fastq.gz -S $filename.sam -p $num_threads -5 $five_trim -3 $three_trim --min-intronlen $min_intl --max-intronlen $max_intl --dta --phred33
            stringtie_non_SRA
            rm $fbname*
            exit
        elif [ $numb -eq 1 ] && [ "$cuff_merge" != 0 ]; then
            filename=$(basename "$f" "_R1_001.fastq.gz")
            hisat2 -x $fbname --rna-strandness $lib_type -1 ${filename}_R1_001.fastq.gz -2 ${filename}_R2_001.fastq.gz -S $filename.sam -p $num_threads -5 $five_trim -3 $three_trim --min-intronlen $min_intl --max-intronlen $max_intl --dta --phred33
            stringtie_non_SRA
            rm $fbname*
            echo "cuffmerge only works with more than one file"
            exit    
        elif [ $numb -gt 1 ] && [ "$cuff_merge" != 0 ]; then
            filename=$(basename "$f" "_R1_001.fastq.gz")
            hisat2 -x $fbname --rna-strandness $lib_type -1 ${filename}_R1_001.fastq.gz -2 ${filename}_R2_001.fastq.gz -S $filename.sam -p $num_threads -5 $five_trim -3 $three_trim --min-intronlen $min_intl --max-intronlen $max_intl --dta --phred33
            stringtie_non_SRA
        fi       
      elif [ "$extension" == "fq.gz" ]; then
        if [ $numb -eq 1 ] && [ "$cuff_merge" == 0 ]; then
            filename=$(basename "$f" "_R1_001.fq.gz")
            hisat2 -x $fbname --rna-strandness $lib_type -1 ${filename}_R1_001.fq.gz -2 ${filename}_R2_001.fq.gz -S $filename.sam -p $num_threads -5 $five_trim -3 $three_trim --min-intronlen $min_intl --max-intronlen $max_intl --dta --phred33
            stringtie_non_SRA
            rm $fbname*
            exit
        elif [ $numb -eq 1 ] && [ "$cuff_merge" != 0 ]; then
            filename=$(basename "$f" "_R1_001.fq.gz")
            hisat2 -x $fbname --rna-strandness $lib_type -1 ${filename}_R1_001.fq.gz -2 ${filename}_R2_001.fq.gz -S $filename.sam -p $num_threads -5 $five_trim -3 $three_trim --min-intronlen $min_intl --max-intronlen $max_intl --dta --phred33
            stringtie_non_SRA
            rm $fbname*
            echo "cuffmerge only works with more than one file"
            exit    
        elif [ $numb -gt 1 ] && [ "$cuff_merge" != 0 ]; then
            filename=$(basename "$f" "_R1_001.fq.gz")
            hisat2 -x $fbname --rna-strandness $lib_type -1 ${filename}_R1_001.fq.gz -2 ${filename}_R2_001.fq.gz -S $filename.sam -p $num_threads -5 $five_trim -3 $three_trim --min-intronlen $min_intl --max-intronlen $max_intl --dta --phred33
            stringtie_non_SRA
        fi       
      elif [ "$extension" != "fastq" ] || [ "$extension" != "fq" ] || [ "$extension" != "fastq.gz" ] || [ "$extension" != "fq.gz" ]; then
        rm $fbname*
        echo "The extension" "$extension" "is not supported. Only .fq, .fq.gz, .fastq, .fastq.gz are only supported" 1>&2        
        exit 64
      fi 
    done
    cuff_merge_fun

# Cufflinks

elif [ ! -z "$left_reads" ] && [ ! -z "$right_reads" ] && [ "$quality_33" != 0 ] && [ "$tra_as" != 0 ] && [ "$tra_cuff" == 0 ]; then
    mkdir "$bam_out"
    numb=$(ls "${left_reads[@]}" | wc -l)
    for f in "${left_reads[@]}"; do
      extension=${f#*.}
      if [ "$extension" == "fastq" ]; then
        if [ $numb -eq 1 ] && [ "$cuff_merge" == 0 ]; then
            filename=$(basename "$f" "_R1_001.fastq")
            hisat2 -x $fbname --rna-strandness $lib_type -1 ${filename}_R1_001.fastq -2 ${filename}_R2_001.fastq -S $filename.sam -p $num_threads -5 $five_trim -3 $three_trim --min-intronlen $min_intl --max-intronlen $max_intl --dta --phred33
            cufflinks_non_SRA
            rm $fbname*
            exit
        elif [ $numb -eq 1 ] && [ "$cuff_merge" != 0 ]; then
            filename=$(basename "$f" "_R1_001.fastq")
            hisat2 -x $fbname --rna-strandness $lib_type -1 ${filename}_R1_001.fastq -2 ${filename}_R2_001.fastq -S $filename.sam -p $num_threads -5 $five_trim -3 $three_trim --min-intronlen $min_intl --max-intronlen $max_intl --dta --phred33
            cufflinks_non_SRA
            rm $fbname*
            echo "cuffmerge only works with more than one file"
            exit    
        elif [ $numb -gt 1 ] && [ "$cuff_merge" != 0 ]; then
            filename=$(basename "$f" "_R1_001.fastq")
            hisat2 -x $fbname --rna-strandness $lib_type -1 ${filename}_R1_001.fastq -2 ${filename}_R2_001.fastq -S $filename.sam -p $num_threads -5 $five_trim -3 $three_trim --min-intronlen $min_intl --max-intronlen $max_intl --dta --phred33
            cufflinks_non_SRA
        fi        
      elif [ "$extension" == "fq" ]; then
        if [ $numb -eq 1 ] && [ "$cuff_merge" == 0 ]; then
            filename=$(basename "$f" "_R1_001.fq")
            hisat2 -x $fbname --rna-strandness $lib_type -1 ${filename}_R1_001.fq -2 ${filename}_R2_001.fq -S $filename.sam -p $num_threads -5 $five_trim -3 $three_trim --min-intronlen $min_intl --max-intronlen $max_intl --dta --phred33
            cufflinks_non_SRA
            rm $fbname*
            exit
        elif [ $numb -eq 1 ] && [ "$cuff_merge" != 0 ]; then
            filename=$(basename "$f" "_R1_001.fq")
            hisat2 -x $fbname --rna-strandness $lib_type -1 ${filename}_R1_001.fq -2 ${filename}_R2_001.fq -S $filename.sam -p $num_threads -5 $five_trim -3 $three_trim --min-intronlen $min_intl --max-intronlen $max_intl --dta --phred33
            cufflinks_non_SRA
            rm $fbname*
            echo "cuffmerge only works with more than one file"
            exit    
        elif [ $numb -gt 1 ] && [ "$cuff_merge" != 0 ]; then
            filename=$(basename "$f" "_R1_001.fq")
            hisat2 -x $fbname --rna-strandness $lib_type -1 ${filename}_R1_001.fq -2 ${filename}_R2_001.fq -S $filename.sam -p $num_threads -5 $five_trim -3 $three_trim --min-intronlen $min_intl --max-intronlen $max_intl --dta --phred33
            cufflinks_non_SRA
        fi       
      elif [ "$extension" == "fastq.gz" ]; then
        if [ $numb -eq 1 ] && [ "$cuff_merge" == 0 ]; then
            filename=$(basename "$f" "_R1_001.fastq.gz")
            hisat2 -x $fbname --rna-strandness $lib_type -1 ${filename}_R1_001.fastq.gz -2 ${filename}_R2_001.fastq.gz -S $filename.sam -p $num_threads -5 $five_trim -3 $three_trim --min-intronlen $min_intl --max-intronlen $max_intl --dta --phred33
            cufflinks_non_SRA
            rm $fbname*
            exit
        elif [ $numb -eq 1 ] && [ "$cuff_merge" != 0 ]; then
            filename=$(basename "$f" "_R1_001.fastq.gz")
            hisat2 -x $fbname --rna-strandness $lib_type -1 ${filename}_R1_001.fastq.gz -2 ${filename}_R2_001.fastq.gz -S $filename.sam -p $num_threads -5 $five_trim -3 $three_trim --min-intronlen $min_intl --max-intronlen $max_intl --dta --phred33
            cufflinks_non_SRA
            rm $fbname*
            echo "cuffmerge only works with more than one file"
            exit    
        elif [ $numb -gt 1 ] && [ "$cuff_merge" != 0 ]; then
            filename=$(basename "$f" "_R1_001.fastq.gz")
            hisat2 -x $fbname --rna-strandness $lib_type -1 ${filename}_R1_001.fastq.gz -2 ${filename}_R2_001.fastq.gz -S $filename.sam -p $num_threads -5 $five_trim -3 $three_trim --min-intronlen $min_intl --max-intronlen $max_intl --dta --phred33
            cufflinks_non_SRA
        fi       
      elif [ "$extension" == "fq.gz" ]; then
        if [ $numb -eq 1 ] && [ "$cuff_merge" == 0 ]; then
            filename=$(basename "$f" "_R1_001.fq.gz")
            hisat2 -x $fbname --rna-strandness $lib_type -1 ${filename}_R1_001.fq.gz -2 ${filename}_R2_001.fq.gz -S $filename.sam -p $num_threads -5 $five_trim -3 $three_trim --min-intronlen $min_intl --max-intronlen $max_intl --dta --phred33
            cufflinks_non_SRA
            rm $fbname*
            exit
        elif [ $numb -eq 1 ] && [ "$cuff_merge" != 0 ]; then
            filename=$(basename "$f" "_R1_001.fq.gz")
            hisat2 -x $fbname --rna-strandness $lib_type -1 ${filename}_R1_001.fq.gz -2 ${filename}_R2_001.fq.gz -S $filename.sam -p $num_threads -5 $five_trim -3 $three_trim --min-intronlen $min_intl --max-intronlen $max_intl --dta --phred33
            cufflinks_non_SRA
            rm $fbname*
            echo "cuffmerge only works with more than one file"
            exit    
        elif [ $numb -gt 1 ] && [ "$cuff_merge" != 0 ]; then
            filename=$(basename "$f" "_R1_001.fq.gz")
            hisat2 -x $fbname --rna-strandness $lib_type -1 ${filename}_R1_001.fq.gz -2 ${filename}_R2_001.fq.gz -S $filename.sam -p $num_threads -5 $five_trim -3 $three_trim --min-intronlen $min_intl --max-intronlen $max_intl --dta --phred33
            cufflinks_non_SRA
        fi       
      elif [ "$extension" != "fastq" ] || [ "$extension" != "fq" ] || [ "$extension" != "fastq.gz" ] || [ "$extension" != "fq.gz" ]; then
        echo "The extension" "$extension" "is not supported. Only .fq, .fq.gz, .fastq, .fastq.gz are only supported" 1>&2
        exit 64
      fi 
    done
    cuff_merge_fun

# Phred 64

elif [ ! -z "$left_reads" ] && [ ! -z "$right_reads" ] && [ "$quality_64" != 0 ] && [ "$tra_as" != 0 ] && [ "$tra_cuff" == 0 ]; then
    mkdir "$bam_out"
    numb=$(ls "${left_reads[@]}" | wc -l)
    for f in "${left_reads[@]}"; do
      extension=${f#*.}
      if [ "$extension" == "fastq" ]; then
        if [ $numb -eq 1 ] && [ "$cuff_merge" == 0 ]; then
            filename=$(basename "$f" "_R1_001.fastq")
            hisat2 -x $fbname --rna-strandness $lib_type -1 ${filename}_R1_001.fastq -2 ${filename}_R2_001.fastq -S $filename.sam -p $num_threads -5 $five_trim -3 $three_trim --min-intronlen $min_intl --max-intronlen $max_intl --dta --phred64
            stringtie_non_SRA
            rm $fbname*
            exit
        elif [ $numb -eq 1 ] && [ "$cuff_merge" != 0 ]; then
            filename=$(basename "$f" "_R1_001.fastq")
            hisat2 -x $fbname --rna-strandness $lib_type -1 ${filename}_R1_001.fastq -2 ${filename}_R2_001.fastq -S $filename.sam -p $num_threads -5 $five_trim -3 $three_trim --min-intronlen $min_intl --max-intronlen $max_intl --dta --phred64
            stringtie_non_SRA
            rm $fbname*
            echo "cuffmerge only works with more than one file"
            exit    
        elif [ $numb -gt 1 ] && [ "$cuff_merge" != 0 ]; then
            filename=$(basename "$f" "_R1_001.fastq")
            hisat2 -x $fbname --rna-strandness $lib_type -1 ${filename}_R1_001.fastq -2 ${filename}_R2_001.fastq -S $filename.sam -p $num_threads -5 $five_trim -3 $three_trim --min-intronlen $min_intl --max-intronlen $max_intl --dta --phred64
            stringtie_non_SRA
        fi        
      elif [ "$extension" == "fq" ]; then
        if [ $numb -eq 1 ] && [ "$cuff_merge" == 0 ]; then
            filename=$(basename "$f" "_R1_001.fq")
            hisat2 -x $fbname --rna-strandness $lib_type -1 ${filename}_R1_001.fq -2 ${filename}_R2_001.fq -S $filename.sam -p $num_threads -5 $five_trim -3 $three_trim --min-intronlen $min_intl --max-intronlen $max_intl --dta --phred64
            stringtie_non_SRA
            rm $fbname*
            exit
        elif [ $numb -eq 1 ] && [ "$cuff_merge" != 0 ]; then
            filename=$(basename "$f" "_R1_001.fq")
            hisat2 -x $fbname --rna-strandness $lib_type -1 ${filename}_R1_001.fq -2 ${filename}_R2_001.fq -S $filename.sam -p $num_threads -5 $five_trim -3 $three_trim --min-intronlen $min_intl --max-intronlen $max_intl --dta --phred64
            stringtie_non_SRA
            rm $fbname*
            echo "cuffmerge only works with more than one file"
            exit    
        elif [ $numb -gt 1 ] && [ "$cuff_merge" != 0 ]; then
            filename=$(basename "$f" "_R1_001.fq")
            hisat2 -x $fbname --rna-strandness $lib_type -1 ${filename}_R1_001.fq -2 ${filename}_R2_001.fq -S $filename.sam -p $num_threads -5 $five_trim -3 $three_trim --min-intronlen $min_intl --max-intronlen $max_intl --dta --phred64
            stringtie_non_SRA
        fi       
      elif [ "$extension" == "fastq.gz" ]; then
        if [ $numb -eq 1 ] && [ "$cuff_merge" == 0 ]; then
            filename=$(basename "$f" "_R1_001.fastq.gz")
            hisat2 -x $fbname --rna-strandness $lib_type -1 ${filename}_R1_001.fastq.gz -2 ${filename}_R2_001.fastq.gz -S $filename.sam -p $num_threads -5 $five_trim -3 $three_trim --min-intronlen $min_intl --max-intronlen $max_intl --dta --phred64
            stringtie_non_SRA
            rm $fbname*
            exit
        elif [ $numb -eq 1 ] && [ "$cuff_merge" != 0 ]; then
            filename=$(basename "$f" "_R1_001.fastq.gz")
            hisat2 -x $fbname --rna-strandness $lib_type -1 ${filename}_R1_001.fastq.gz -2 ${filename}_R2_001.fastq.gz -S $filename.sam -p $num_threads -5 $five_trim -3 $three_trim --min-intronlen $min_intl --max-intronlen $max_intl --dta --phred64
            stringtie_non_SRA
            rm $fbname*
            echo "cuffmerge only works with more than one file"
            exit    
        elif [ $numb -gt 1 ] && [ "$cuff_merge" != 0 ]; then
            filename=$(basename "$f" "_R1_001.fastq.gz")
            hisat2 -x $fbname --rna-strandness $lib_type -1 ${filename}_R1_001.fastq.gz -2 ${filename}_R2_001.fastq.gz -S $filename.sam -p $num_threads -5 $five_trim -3 $three_trim --min-intronlen $min_intl --max-intronlen $max_intl --dta --phred64
            stringtie_non_SRA
        fi       
      elif [ "$extension" == "fq.gz" ]; then
        if [ $numb -eq 1 ] && [ "$cuff_merge" == 0 ]; then
            filename=$(basename "$f" "_R1_001.fq.gz")
            hisat2 -x $fbname --rna-strandness $lib_type -1 ${filename}_R1_001.fq.gz -2 ${filename}_R2_001.fq.gz -S $filename.sam -p $num_threads -5 $five_trim -3 $three_trim --min-intronlen $min_intl --max-intronlen $max_intl --dta --phred64
            stringtie_non_SRA
            rm $fbname*
            exit
        elif [ $numb -eq 1 ] && [ "$cuff_merge" != 0 ]; then
            filename=$(basename "$f" "_R1_001.fq.gz")
            hisat2 -x $fbname --rna-strandness $lib_type -1 ${filename}_R1_001.fq.gz -2 ${filename}_R2_001.fq.gz -S $filename.sam -p $num_threads -5 $five_trim -3 $three_trim --min-intronlen $min_intl --max-intronlen $max_intl --dta --phred64
            stringtie_non_SRA
            rm $fbname*
            echo "cuffmerge only works with more than one file"
            exit    
        elif [ $numb -gt 1 ] && [ "$cuff_merge" != 0 ]; then
            filename=$(basename "$f" "_R1_001.fq.gz")
            hisat2 -x $fbname --rna-strandness $lib_type -1 ${filename}_R1_001.fq.gz -2 ${filename}_R2_001.fq.gz -S $filename.sam -p $num_threads -5 $five_trim -3 $three_trim --min-intronlen $min_intl --max-intronlen $max_intl --dta --phred64
            stringtie_non_SRA
        fi       
      elif [ "$extension" != "fastq" ] || [ "$extension" != "fq" ] || [ "$extension" != "fastq.gz" ] || [ "$extension" != "fq.gz" ]; then
        echo "The extension" "$extension" "is not supported. Only .fq, .fq.gz, .fastq, .fastq.gz are only supported" 1>&2
        exit 64
      fi 
    done
    cuff_merge_fun

elif [ ! -z "$left_reads" ] && [ ! -z "$right_reads" ] && [ "$quality_64" != 0 ] && [ "$tra_as" != 0 ] && [ "$tra_cuff" == 0 ]; then
    mkdir "$bam_out"
    numb=$(ls "${left_reads[@]}" | wc -l)
    for f in "${left_reads[@]}"; do
      extension=${f#*.}
      if [ "$extension" == "fastq" ]; then
        if [ $numb -eq 1 ] && [ "$cuff_merge" == 0 ]; then
            filename=$(basename "$f" "_R1_001.fastq")
            hisat2 -x $fbname --rna-strandness $lib_type -1 ${filename}_R1_001.fastq -2 ${filename}_R2_001.fastq -S $filename.sam -p $num_threads -5 $five_trim -3 $three_trim --min-intronlen $min_intl --max-intronlen $max_intl --dta --phred64
            cufflinks_non_SRA
            rm $fbname*
            exit
        elif [ $numb -eq 1 ] && [ "$cuff_merge" != 0 ]; then
            filename=$(basename "$f" "_R1_001.fastq")
            hisat2 -x $fbname --rna-strandness $lib_type -1 ${filename}_R1_001.fastq -2 ${filename}_R2_001.fastq -S $filename.sam -p $num_threads -5 $five_trim -3 $three_trim --min-intronlen $min_intl --max-intronlen $max_intl --dta --phred64
            cufflinks_non_SRA
            rm $fbname*
            echo "cuffmerge only works with more than one file"
            exit    
        elif [ $numb -gt 1 ] && [ "$cuff_merge" != 0 ]; then
            filename=$(basename "$f" "_R1_001.fastq")
            hisat2 -x $fbname --rna-strandness $lib_type -1 ${filename}_R1_001.fastq -2 ${filename}_R2_001.fastq -S $filename.sam -p $num_threads -5 $five_trim -3 $three_trim --min-intronlen $min_intl --max-intronlen $max_intl --dta --phred64
            cufflinks_non_SRA
        fi        
      elif [ "$extension" == "fq" ]; then
        if [ $numb -eq 1 ] && [ "$cuff_merge" == 0 ]; then
            filename=$(basename "$f" "_R1_001.fq")
            hisat2 -x $fbname --rna-strandness $lib_type -1 ${filename}_R1_001.fq -2 ${filename}_R2_001.fq -S $filename.sam -p $num_threads -5 $five_trim -3 $three_trim --min-intronlen $min_intl --max-intronlen $max_intl --dta --phred64
            cufflinks_non_SRA
            rm $fbname*
            exit
        elif [ $numb -eq 1 ] && [ "$cuff_merge" != 0 ]; then
            filename=$(basename "$f" "_R1_001.fq")
            hisat2 -x $fbname --rna-strandness $lib_type -1 ${filename}_R1_001.fq -2 ${filename}_R2_001.fq -S $filename.sam -p $num_threads -5 $five_trim -3 $three_trim --min-intronlen $min_intl --max-intronlen $max_intl --dta --phred64
            cufflinks_non_SRA
            rm $fbname*
            echo "cuffmerge only works with more than one file"
            exit    
        elif [ $numb -gt 1 ] && [ "$cuff_merge" != 0 ]; then
            filename=$(basename "$f" "_R1_001.fq")
            hisat2 -x $fbname --rna-strandness $lib_type -1 ${filename}_R1_001.fq -2 ${filename}_R2_001.fq -S $filename.sam -p $num_threads -5 $five_trim -3 $three_trim --min-intronlen $min_intl --max-intronlen $max_intl --dta --phred64
            cufflinks_non_SRA
        fi       
      elif [ "$extension" == "fastq.gz" ]; then
        if [ $numb -eq 1 ] && [ "$cuff_merge" == 0 ]; then
            filename=$(basename "$f" "_R1_001.fastq.gz")
            hisat2 -x $fbname --rna-strandness $lib_type -1 ${filename}_R1_001.fastq.gz -2 ${filename}_R2_001.fastq.gz -S $filename.sam -p $num_threads -5 $five_trim -3 $three_trim --min-intronlen $min_intl --max-intronlen $max_intl --dta --phred64
            cufflinks_non_SRA
            rm $fbname*
            exit
        elif [ $numb -eq 1 ] && [ "$cuff_merge" != 0 ]; then
            filename=$(basename "$f" "_R1_001.fastq.gz")
            hisat2 -x $fbname --rna-strandness $lib_type -1 ${filename}_R1_001.fastq.gz -2 ${filename}_R2_001.fastq.gz -S $filename.sam -p $num_threads -5 $five_trim -3 $three_trim --min-intronlen $min_intl --max-intronlen $max_intl --dta --phred64
            cufflinks_non_SRA
            rm $fbname*
            echo "cuffmerge only works with more than one file"
            exit    
        elif [ $numb -gt 1 ] && [ "$cuff_merge" != 0 ]; then
            filename=$(basename "$f" "_R1_001.fastq.gz")
            hisat2 -x $fbname --rna-strandness $lib_type -1 ${filename}_R1_001.fastq.gz -2 ${filename}_R2_001.fastq.gz -S $filename.sam -p $num_threads -5 $five_trim -3 $three_trim --min-intronlen $min_intl --max-intronlen $max_intl --dta --phred64
            cufflinks_non_SRA
        fi       
      elif [ "$extension" == "fq.gz" ]; then
        if [ $numb -eq 1 ] && [ "$cuff_merge" == 0 ]; then
            filename=$(basename "$f" "_R1_001.fq.gz")
            hisat2 -x $fbname --rna-strandness $lib_type -1 ${filename}_R1_001.fq.gz -2 ${filename}_R2_001.fq.gz -S $filename.sam -p $num_threads -5 $five_trim -3 $three_trim --min-intronlen $min_intl --max-intronlen $max_intl --dta --phred64
            cufflinks_non_SRA
            rm $fbname*
            exit
        elif [ $numb -eq 1 ] && [ "$cuff_merge" != 0 ]; then
            filename=$(basename "$f" "_R1_001.fq.gz")
            hisat2 -x $fbname --rna-strandness $lib_type -1 ${filename}_R1_001.fq.gz -2 ${filename}_R2_001.fq.gz -S $filename.sam -p $num_threads -5 $five_trim -3 $three_trim --min-intronlen $min_intl --max-intronlen $max_intl --dta --phred64
            cufflinks_non_SRA
            rm $fbname*
            echo "cuffmerge only works with more than one file"
            exit    
        elif [ $numb -gt 1 ] && [ "$cuff_merge" != 0 ]; then
            filename=$(basename "$f" "_R1_001.fq.gz")
            hisat2 -x $fbname --rna-strandness $lib_type -1 ${filename}_R1_001.fq.gz -2 ${filename}_R2_001.fq.gz -S $filename.sam -p $num_threads -5 $five_trim -3 $three_trim --min-intronlen $min_intl --max-intronlen $max_intl --dta --phred64
            cufflinks_non_SRA
        fi       
      elif [ "$extension" != "fastq" ] || [ "$extension" != "fq" ] || [ "$extension" != "fastq.gz" ] || [ "$extension" != "fq.gz" ]; then
        echo "The extension" "$extension" "is not supported. Only .fq, .fq.gz, .fastq, .fastq.gz are only supported" 1>&2
        exit 64
      fi 
    done
    cuff_merge_fun

# ############################################################################################################################################################################################################################
# # Single end reads
# ############################################################################################################################################################################################################################

# Phred 33

# Stringtie
elif [ ! -z "$single_reads" ] && [ "$quality_33" != 0 ] && [ "$tra_as" != 0 ] && [ "$tra_cuff" == 0 ]; then
    mkdir "$bam_out"
    numb=$(ls "${single_reads[@]}" | wc -l)
    for f in "${single_reads[@]}"; do
      if [ $numb -eq 1 ] && [ "$cuff_merge" == 0 ]; then
        filename=$(basename "$f")
        hisat2 -x $fbname --rna-strandness $lib_type -U $f -S $filename.sam -p $num_threads -5 $five_trim -3 $three_trim --min-intronlen $min_intl --max-intronlen $max_intl --dta --phred33
        stringtie_non_SRA
        rm $fbname*
        exit
      elif [ $numb -eq 1 ] && [ "$cuff_merge" != 0 ]; then
        filename=$(basename "$f")
        hisat2 -x $fbname --rna-strandness $lib_type -U $f -S $filename.sam -p $num_threads -5 $five_trim -3 $three_trim --min-intronlen $min_intl --max-intronlen $max_intl --dta --phred33
        stringtie_non_SRA
        rm $fbname*
        echo "cuffmerge only works with more than one file"
        exit   
      elif [ $numb -gt 1 ] && [ "$cuff_merge" != 0 ]; then
        filename=$(basename "$f")
        hisat2 -x $fbname --rna-strandness $lib_type -U $f -S $filename.sam -p $num_threads -5 $five_trim -3 $three_trim --min-intronlen $min_intl --max-intronlen $max_intl --dta --phred33
        stringtie_non_SRA
      fi 
    done
    cuff_merge_fun 

# Cufflinks
elif [ ! -z "$single_reads" ] && [ "$quality_33" != 0 ] && [ "$tra_as" == 0 ] && [ "$tra_cuff" != 0 ]; then
    mkdir "$bam_out"
    numb=$(ls "${single_reads[@]}" | wc -l)
    for f in "${single_reads[@]}"; do
      if [ $numb -eq 1 ] && [ "$cuff_merge" == 0 ]; then 
        filename=$(basename "$f")
        hisat2 -x $fbname --rna-strandness $lib_type -U $f -S $filename.sam -p $num_threads -5 $five_trim -3 $three_trim --min-intronlen $min_intl --max-intronlen $max_intl --dta-cufflinks --phred33
        cufflinks_non_SRA
        rm $fbname*
        exit
      elif [ $numb -eq 1 ] && [ "$cuff_merge" != 0 ]; then 
        filename=$(basename "$f")
        hisat2 -x $fbname --rna-strandness $lib_type -U $f -S $filename.sam -p $num_threads -5 $five_trim -3 $three_trim --min-intronlen $min_intl --max-intronlen $max_intl --dta-cufflinks --phred33
        cufflinks_non_SRA
        rm $fbname*
        echo "cuffmerge only works with more than one file"
        exit  
      elif [ $numb -gt 1 ] && [ "$cuff_merge" != 0 ]; then
        filename=$(basename "$f")
        hisat2 -x $fbname --rna-strandness $lib_type -U $f -S $filename.sam -p $num_threads -5 $five_trim -3 $three_trim --min-intronlen $min_intl --max-intronlen $max_intl --dta-cufflinks --phred33
        cufflinks_non_SRA
      fi          
    done
    cuff_merge_fun

# Phred 64

# Stringtie
elif [ ! -z "$single_reads" ] && [ "$quality_64" != 0 ] && [ "$tra_as" != 0 ] && [ "$tra_cuff" == 0 ]; then
    mkdir "$bam_out"
    numb=$(ls "${single_reads[@]}" | wc -l)
    for f in "${single_reads[@]}"; do
      if [ $numb -eq 1 ] && [ "$cuff_merge" == 0 ]; then
        filename=$(basename "$f")
        hisat2 -x $fbname --rna-strandness $lib_type -U $f -S $filename.sam -p $num_threads -5 $five_trim -3 $three_trim --min-intronlen $min_intl --max-intronlen $max_intl --dta --phred64
        stringtie_non_SRA
        rm $fbname*
        exit 
      elif [ $numb -eq 1 ] && [ "$cuff_merge" != 0 ]; then
        filename=$(basename "$f")
        hisat2 -x $fbname --rna-strandness $lib_type -U $f -S $filename.sam -p $num_threads -5 $five_trim -3 $three_trim --min-intronlen $min_intl --max-intronlen $max_intl --dta --phred64
        stringtie_non_SRA
        rm $fbname*
        echo "cuffmerge only works with more than one file"
        exit  
      elif [ $numb -gt 1 ] && [ "$cuff_merge" != 0 ]; then
        filename=$(basename "$f")
        hisat2 -x $fbname --rna-strandness $lib_type -U $f -S $filename.sam -p $num_threads -5 $five_trim -3 $three_trim --min-intronlen $min_intl --max-intronlen $max_intl --dta --phred64
        stringtie_non_SRA
      fi 
    done
    cuff_merge_fun 

# Cufflinks
elif [ ! -z "$single_reads" ] && [ "$quality_64" != 0 ] && [ "$tra_as" == 0 ] && [ "$tra_cuff" != 0 ]; then
    mkdir "$bam_out"
    numb=$(ls "${single_reads[@]}" | wc -l)
    for f in "${single_reads[@]}"; do
      if [ $numb -eq 1 ] && [ "$cuff_merge" == 0 ]; then 
        filename=$(basename "$f")
        hisat2 -x $fbname --rna-strandness $lib_type -U $f -S $filename.sam -p $num_threads -5 $five_trim -3 $three_trim --min-intronlen $min_intl --max-intronlen $max_intl --dta-cufflinks --phred64
        cufflinks_non_SRA
        rm $fbname*
        exit
      elif [ $numb -eq 1 ] && [ "$cuff_merge" != 0 ]; then 
        filename=$(basename "$f")
        hisat2 -x $fbname --rna-strandness $lib_type -U $f -S $filename.sam -p $num_threads -5 $five_trim -3 $three_trim --min-intronlen $min_intl --max-intronlen $max_intl --dta-cufflinks --phred64
        cufflinks_non_SRA
        rm $fbname*
        echo "cuffmerge only works with more than one file"
        exit  
      elif [ $numb -gt 1 ] && [ "$cuff_merge" != 0 ]; then
        filename=$(basename "$f")
        hisat2 -x $fbname --rna-strandness $lib_type -U $f -S $filename.sam -p $num_threads -5 $five_trim -3 $three_trim --min-intronlen $min_intl --max-intronlen $max_intl --dta-cufflinks --phred64
        cufflinks_non_SRA
      fi          
    done
    cuff_merge_fun

# ############################################################################################################################################################################################################################
# # SRA reads
# ############################################################################################################################################################################################################################
# phred 33

# Stringtie
elif [ ! -z $sra_id ] && [ "$quality_33" != 0 ] && [ "$tra_as" != 0 ] && [ "$tra_cuff" == 0 ]; then
  if [[ -f $sra_id ]]; then
    mkdir "$bam_out"    
    while read f; do
      hisat2 -x $fbname --rna-strandness $lib_type --sra-acc $f -S $f.sam -p $num_threads -5 $five_trim -3 $three_trim --min-intronlen $min_intl --max-intronlen $max_intl --dta --phred33
      stringtie_SRA_multi      
      done < "$sra_id" 
      cuff_merge_fun
  else    
    mkdir "$bam_out"
    hisat2 -x $fbname --rna-strandness $lib_type --sra-acc $sra_id -S $sra_id.sam -p $num_threads -5 $five_trim -3 $three_trim --min-intronlen $min_intl --max-intronlen $max_intl --dta --phred33
    stringtie_SRA_single
    if [ "$cuff_merge" != 0 ]; then
      echo "cuffmerge only works with more than one SRA accesions. Use File containing SRA id's option to user cuffmerge"
    fi  
  fi

# Cufflinks
elif [ ! -z $sra_id ] && [ "$quality_33" != 0 ] && [ "$tra_as" == 0 ] && [ "$tra_cuff" != 0 ]; then
  if [[ -f $sra_id ]]; then
    mkdir "$bam_out"    
    while read f; do
      hisat2 -x $fbname --rna-strandness $lib_type --sra-acc $f -S $f.sam -p $num_threads -5 $five_trim -3 $three_trim --min-intronlen $min_intl --max-intronlen $max_intl --dta-cufflinks --phred33
      cufflinks_SRA_multi      
      done < "$sra_id"
      cuff_merge_fun
  else    
    mkdir "$bam_out"
    hisat2 -x $fbname --rna-strandness $lib_type --sra-acc $sra_id -S $sra_id.sam -p $num_threads -5 $five_trim -3 $three_trim --min-intronlen $min_intl --max-intronlen $max_intl --dta-cufflinks --phred33
    cufflinks_SRA_single
    if [ "$cuff_merge" != 0 ]; then
      echo "cuffmerge only works with more than one SRA accesions. Use File containing SRA id's option to user cuffmerge"
    fi
  fi

# Phred 64

# Stringtie
elif [ ! -z $sra_id ] && [ "$quality_64" != 0 ] && [ "$tra_as" != 0 ] && [ "$tra_cuff" == 0 ]; then
  if [[ -f $sra_id ]]; then
    mkdir "$bam_out"    
    while read f; do
      hisat2 -x $fbname --rna-strandness $lib_type --sra-acc $f -S $f.sam -p $num_threads -5 $five_trim -3 $three_trim --min-intronlen $min_intl --max-intronlen $max_intl --dta --phred64
      stringtie_SRA_multi      
      done < "$sra_id"
      cuff_merge_fun
  else    
    mkdir "$bam_out"
    hisat2 -x $fbname --rna-strandness $lib_type --sra-acc $sra_id -S $sra_id.sam -p $num_threads -5 $five_trim -3 $three_trim --min-intronlen $min_intl --max-intronlen $max_intl --dta --phred64
    stringtie_SRA_single
    if [ "$cuff_merge" != 0 ]; then
      echo "cuffmerge only works with more than one SRA accesions. Use File containing SRA id's option to user cuffmerge"
    fi
  fi 

# Cufflinks
elif [ ! -z $sra_id ] && [ "$quality_64" != 0 ] && [ "$tra_as" == 0 ] && [ "$tra_cuff" != 0 ]; then
  if [[ -f $sra_id ]]; then
    mkdir "$bam_out"    
    while read f; do
      hisat2 -x $fbname --rna-strandness $lib_type --sra-acc $f -S $f.sam -p $num_threads -5 $five_trim -3 $three_trim --min-intronlen $min_intl --max-intronlen $max_intl --dta-cufflinks --phred64
      cufflinks_SRA_multi      
      done < "$sra_id"
      cuff_merge_fun
  else    
    mkdir "$bam_out"
    hisat2 -x $fbname --rna-strandness $lib_type --sra-acc $sra_id -S $sra_id.sam -p $num_threads -5 $five_trim -3 $three_trim --min-intronlen $min_intl --max-intronlen $max_intl --dta-cufflinks --phred64
    cufflinks_SRA_single
    if [ "$cuff_merge" != 0 ]; then
      echo "cuffmerge only works with more than one SRA accesions. Use File containing SRA id's option to user cuffmerge"
    fi
  fi
fi

# Clean up the reference genomes
rm $fbname*
