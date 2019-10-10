#!/bin/sh

## command line arguments:
##  -i	input_file
##  -o	output_file
##	-n	number of lines
##	-r	invert (extract all but first n lines)

INPUT=
OUTPUT=tail_output.txt
LINES=5
FLAGS="-n"

while getopts “ri:o:n:” OPTION
do
     case $OPTION in
		i)
			INPUT=$OPTARG
			;;
		o)
			OUTPUT=$OPTARG
			;;
		n)
			LINES=$OPTARG
			;;
		r)
			FLAGS="${FLAGS}+"
			;;
	esac
done

set -x
tail "${FLAGS}${LINES}" $INPUT > $OUTPUT
set -x
