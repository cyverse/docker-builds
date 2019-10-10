#!/bin/bash
# Usage: evalTicket.sh [ticket] [target] [command]
# ticket: iRODS hex code verification
# target: directory for extracting/exporting data
# command: 'iget' to extract data (iget) 'iput' to export data (iput)

# Parse iRODS input file to get ticket and data
var=$2
var=$(echo ${var:0:1})
type=$1

case $var in
'#')
        ;;
*)
    case $type in
    iput)
        ticket=$(echo $2 | cut -f 1 -d ,)
        target=$(echo $2 | cut -f 2 -d ,)
        username=$(echo $target | cut -f 4 -d /)
        sourcedir=$PWD/loadingdock/userdata
        sourcefldr=dataout
        printf "\n\nUSER: ${username}\n\n"
        printf "Using ticket ${ticket} to ${type} ${source} into ${target}\n"
        #iput -v -t $ticket -r $source $target

        # Upload output structure to CyVerse
        # 1) Create master output directory into user's CyVerse directory
        # 2) Grant ownership permission to the `job` user in case inheritance is enabled.
        # 3) List master directory structure then call xargs to recursively iput each
        #       file/folder into user's master directory
        # 4) Grant user permissions to master directory
        # 5) Revoke anonymous user's permissions to directory
        iput -Vrt "$ticket" "$sourcedir/$sourcefldr" "$target" > /dev/null 2>&1
        ichmod -r own job "$target/$sourcefldr"
        for file in "$sourcedir/$sourcefldr"/*; do
            iput -Vr "$file" "$target/$sourcefldr/"
        done

        ichmod -r own $username  "$target/$sourcefldr"
        hostname=job
        ichmod -r null $hostname "$target/$sourcefldr"
        ;;

    iget)
        ticket=$(echo $2 | cut -f 1 -d ,)
        source=$(echo $2 | cut -f 2 -d ,)
        target=$PWD/loadingdock/userdata/datain/
        printf "Using ticket ${ticket} to ${type} ${source} into ${target}\n"
        iget -v -t "$ticket" -r "$source $target"
        ;;

    *)
        ;;
    esac
esac
