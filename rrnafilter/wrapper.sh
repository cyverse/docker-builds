#!/usr/bin/bash

while getopts i:r: option
do
        case "${option}"
        in

                i) input=${OPTARG};;
		r) filter=${OPTARG};;
        esac
done

ARGS=''

#IF STATEMENTS EXIST FOR EACH OPTIONAL PARAMETER
if [ -n "${filter}" ]; then ARGS="$ARGS -r $filter"; fi

cd /usr/local/bin/

java -Xmx4g -jar /usr/local/bin/rRNAFilter_commandline.jar -i /de-app-work/$input $ARGS


