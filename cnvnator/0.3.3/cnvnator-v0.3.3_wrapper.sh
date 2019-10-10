#!/bin/bash
# Upendra Kumar Devisetty
# Wrapper script to identify CNV from bam file

usage() {
      echo ""
      echo "Usage : sh $0 -c chromosome -u -b bamfile(s) -r referencegenome -i histogram_bin_size -s stat_bin_size -p partition_bin_size -a call_bin_size -f prefix -o output"
      echo ""

cat <<'EOF'

  -c <chromosome id>

  -u unique to get the correct fraction of reads mapped with q0 quality

  -b </path/to/bam file(s)>

  -r </path/to/reference genome file>

  -i <histogram bin size>

  -s <stat bin size>

  -p <partition bin size>

  -a <call bin size>

  -f <prefix>
  
  -o <output>

  -h Show this usage information

EOF
    exit 0
}

unique=0

while getopts ":c:hb:ur:i:s:p:a:f:o:" opt; do
  case $opt in
    c)
     chromosome=$OPTARG
      ;;
    u)
     unique=$OPTARG
      ;;
    b)
     bamfile+=("$OPTARG")
      ;;
    h)
     usage
     exit 1
      ;;    
    r)
     referencegenome=$OPTARG
      ;;
    i)
     hist=$OPTARG
      ;;
    s)
     stats=$OPTARG
      ;;
    p)
     part=$OPTARG
      ;;  
    a)
     call=$OPTARG
      ;;
    f)
     prefix=$OPTARG
      ;;
    o)
     output=$OPTARG
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

if [[ -f $chromosome ]]; then
   
  new=$(tr '\n' " " < $chromosome)

  cat $chromosome | while read line

  do 
    echo $line > $line.txt

    perl /subset-fasta.pl $line.txt $referencegenome

    mv $line.txt.fa $line.fa

  done
  
  if [ "$unique" != 0 ]; then

    for val in "${bamfile[@]}"

      do

        cnvnator -root cnvnator.root -unique -chrom $new -tree $val 

      done

  fi

  cnvnator -root cnvnator.root -genome $referencegenome -chrom $new -his $hist

  cnvnator -root cnvnator.root -genome $referencegenome -chrom $new -stat $stats

  cnvnator -root cnvnator.root -genome $referencegenome -chrom $new -partition $part
  
  cnvnator -root cnvnator.root -genome $referencegenome -chrom $new -call $call > $output."cnvnator"

else

    echo $chromosome > $chromosome.txt

    perl /subset-fasta.pl $chromosome.txt $referencegenome

    mv $chromosome.txt.fa $chromosome.fa

  if [ "$unique" != 0 ]; then

    for val in "${bamfile[@]}"

      do

        cnvnator -root cnvnator.root -unique -chrom $chromosome -tree $val 

      done

  fi

  cnvnator -root cnvnator.root -genome $referencegenome -chrom $chromosome -his $hist

  cnvnator -root cnvnator.root -genome $referencegenome -chrom $chromosome -stat $stats

  cnvnator -root cnvnator.root -genome $referencegenome -chrom $chromosome -partition $part    

  cnvnator -root cnvnator.root -genome $referencegenome -chrom $chromosome -call $call > $output."cnvnator"

fi

perl /cnvnator2VCF.pl -prefix $prefix $output."cnvnator" > $output.vcf 

rm *.txt *.fa
