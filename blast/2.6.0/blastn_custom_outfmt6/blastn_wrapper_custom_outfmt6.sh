#!/bin/bash

set +o xtrace

set -e
set -x

database="$1"

INPUT_RD=$(basename $database)

if [ -r $INPUT_RD/*.nal ] ; then 
   dbname=$(basename $INPUT_RD/*.nsq .nsq);
else
   dbname=$(basename $INPUT_RD/*.nhr .nhr)
fi

if [ -z "$INPUT_RD" ]; then
echo "Database $database is improperly formatted";
exit 1;
fi

cp $INPUT_RD/$dbname* ./

query="$2"
outfilename="$3"
outformat="$4"
evalue="$5"
options1="$7"
options2="$8"

if [ "$6" == ungapped ]; then options="-ungapped"; fi

if [ "$4" == title ]; then 
	outformat="6 std sscinames stitle"
	wget ftp://ftp.ncbi.nlm.nih.gov/blast/db/taxdb.tar.gz
	tar xvf taxdb.tar.gz
fi

blastn -num_threads 4 -db $dbname -query $query -out $outfilename -outfmt "$outformat" -evalue $evalue $options $options1 $options2

rm -r $dbname*.n*

rm taxdb*

