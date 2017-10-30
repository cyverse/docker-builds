#!/bin/bash

cat $3 | sed 's/ILLUMINACLIP:/ILLUMINACLIP:\/staging\//' > trimfile.txt
mkdir trimout
cuts=`cat trimfile.txt | grep -v \#`

if [ $1 == PE ]; then
        for x in $2/*
        do

        if [[ "$x" =~ .*_R1_001\.fq.gz$ ]]; 
             then
             z=$(basename $x _R1_001.fq.gz)
             java -jar -Xmx1024m /root/Trimmomatic-0.33/trimmomatic-0.33.jar PE "$2"/"$z"'_R1_001.fq.gz' "$2"/"$z"'_R2_001.fq.gz' trimout/'trmPr_'"$z"'_R1_001.fq.gz' trimout/'trmS_'"$z"'_R1_001.fq.gz' trimout/'trmPr_'"$z"'_R2_001.fq.gz' trimout/'trmS_'"$z"'_R2_001.fq.gz' ILLUMINACLIP:$4:2:30:10 $cuts
	     mv trimfile.txt trimmersettings.txt
	     exit
               
        elif [[ "$x" =~ .*_R1_001\.fastq.gz$ ]]; 
             then
             z=$(basename $x _R1_001.fastq.gz)
             java -jar -Xmx1024m /root/Trimmomatic-0.33/trimmomatic-0.33.jar PE "$2"/"$z"'_R1_001.fastq.gz' "$2"/"$z"'_R2_001.fastq.gz' trimout/'trmPr_'"$z"'_R1_001.fastq.gz' trimout/'trmS_'"$z"'_R1_001.fastq.gz' trimout/'trmPr_'"$z"'_R2_001.fastq.gz' trimout/'trmS_'"$z"'_R2_001.fastq.gz' ILLUMINACLIP:$4:2:30:10 $cuts
	     mv trimfile.txt trimmersettings.txt
             exit

	else
	     rm -r trimout trimfile.txt

             echo "The input folder - $2 should contain reads that have extension _R1_001.fq.gz/_R2_001.fq.gz or _R1_001.fastq.gz/_R1_001.fastq.gz

		   Example: EP1_C1_CGATGT_L004_R1_001.fq.gz and EP1_C1_CGATGT_L004_R2_001.fq.gz

		   Please rename your reads and try again" 1>&2

             exit 64

	fi
	done

elif [ $1 == SE ]; then
        for x in $2/*
                do
        	y=$(basename $x)
        	z='trimS'$y
                java -jar -Xmx1024m /root/Trimmomatic-0.33/trimmomatic-0.33.jar SE $x trimout/$z ILLUMINACLIP:$4:2:30:10 $cuts
		mv trimfile.txt trimmersettings.txt
                done
fi
