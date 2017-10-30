#!/bin/bash

set +o xtrace
set -u
set -e

database="$1"

INPUT_RD=$(basename $database)

dbtype="$2"

if [ $dbtype = "nucl" ]; then
	if [ -r $INPUT_RD/*.nin ] ; then 
   	   dbname=$(basename $INPUT_RD/*.nin .nin);
        else
           dbname=$(basename $INPUT_RD/*.nsq .nsq)
	fi
elif [ $dbtype = "prot" ]; then
        if [ -r $INPUT_RD/*.pin ] ; then
           dbname=$(basename $INPUT_RD/*.pin .pin);
        else
           dbname=$(basename $INPUT_RD/*.psq .psq)
	fi
fi

if [ -z "$INPUT_RD" ]; then
echo "Database $database is improperly formatted";
exit 1;
fi

cp $INPUT_RD/$dbname* ./

entry="$3"
out="$4"
file_type="$5"

if [ $file_type == "batch" ]; then

 blastdbcmd -db $dbname -entry_batch $entry -out $out -dbtype $dbtype

elif [ $file_type == "single" ]; then

 range="$6"
 blastdbcmd -db $dbname -entry $entry -out $out -dbtype $dbtype -range $range

fi

rm blastdb_p.p*
