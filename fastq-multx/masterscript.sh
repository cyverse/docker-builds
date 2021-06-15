#!/bin/bash

#######################################################################################################
##SET UP OPTIONS

while getopts B:bd:eG:g:HI:i:L:l:m:no:q:R:r:v:x option
do
        case "${option}"
        in

		B) fromfile=${OPTARG};;
                b) forcebegin=true;;
		d) mindist=${OPTARG};;
                e) forceend=true;;
		G) group=${OPTARG};;
                g) freq=${OPTARG};;
		H) header=true;;
		i) indexinputs+=${OPTARG};;
		m) mismatch=${OPTARG};;
		n) codes=true;;
		o) out=${OPTARG};;
                l) master=${OPTARG};;
		L) readlist=${OPTARG};;
		q) minphred=${OPTARG};;
		R) read1inputs+=${OPTARG};;
		r) read2inputs+=${OPTARG};;
		t) threshold=${OPTARG};;
		v) mate=${OPTARG};;
		x) trim=true;;
        esac
done
#####################################################################################################

ARGS=''
OARGS=''


#IF STATEMENTS EXIST FOR EACH OPTIONAL PARAMETER


if [ -n "${threshold}" ]; then ARGS="$ARGS -t $threshold"; fi
if [ -n "${minphred}" ]; then ARGS="$ARGS -q $minphred"; fi
if [ -n "${readlist}" ]; then ARGS="$ARGS -L $readlist"; fi
if [ -n "${group}" ]; then ARGS="$ARGS -G $group"; fi
if [ -n "${mindist}" ]; then ARGS="$ARGS -d $mindist"; fi
if [ -n "${fromfile}" ]; then ARGS="$ARGS -B $fromfile"; fi
if [ -n "${freq}" ]; then ARGS="$ARGS -g $freq"; fi
if [ -n "${master}" ]; then ARGS="$ARGS -l $master"; fi
if [ -n "${mismatch}" ]; then ARGS="$ARGS -m $mismatch"; fi
if [ -n "${mate}" ]; then ARGS="$ARGS -v $mate"; fi
if [[ "$forcebegin" = "true" ]]; then ARGS="$ARGS -b"; fi
if [[ "$forceend" = "true" ]]; then ARGS="$ARGS -e"; fi
if [[ "$trim" = "true" ]]; then ARGS="$ARGS -x"; fi
if [[ "$codes" = "true" ]]; then ARGS="$ARGS -n"; fi
if [[ "$header" = "true" ]]; then ARGS="$ARGS -H"; fi

for k in "${indexinputs[@]}"; do
	ARGS="$ARGS $k -o n/a"
done


for i in "${read1inputs[@]}"; do
	basename1=$(basename $i | cut -d. -f1)
	ARGS="$ARGS $i" 
	OARGS="$OARGS -o '$basename1'_%_R1.fastq" 
done


for j in "${read2inputs[@]}"; do
	basename2=$(basename $j | cut -d. -f1)
	ARGS="$ARGS $j" 
	OARGS="$OARGS -o '$basename2'_%_R2.fastq" 
done


fastq-multx $ARGS $OARGS
