#!/bin/bash

set +o xtrace

database="$1"

INPUT_RD=$(basename $database)

if [ -r $INPUT_RD/*.pal ] ; then 
	dbname=$(basename $INPUT_RD/*.pal .pal);
else
	dbname=$(basename $INPUT_RD/*.psq .psq)
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
options1="$6"
options2="$7"

blastp -num_threads 4 -db $dbname -query $query -out $outfilename -outfmt $outformat -evalue $evalue $options1 $options2

rm $dbname*.p*
