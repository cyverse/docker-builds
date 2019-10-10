#!/bin/bash
fin=$1
icmd=$2

# Determine if extracting or exporting data
case $icmd in
    iget)

        source=$PWD/loadingdock/userdata/datain/
        ;;

    iput)
        source=$PWD/loadingdock/userdata/dataout/
        ;;

    *)
        echo "error with icommand, ..."
        source=$PWD/loadingdock/userdata/dataout/
        ;;

esac

# Parse through ticket list to either extract or export data
while IFS='' read -r line || [[ -n "$line" ]]; do
    bash evalTicket.sh $icmd $line
done < "$fin"
