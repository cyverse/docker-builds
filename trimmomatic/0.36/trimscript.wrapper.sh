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
             java -jar -Xmx1024m /root/Trimmomatic-0.36/trimmomatic-0.36.jar PE "$2"/"$z"'_R1_001.fq.gz' "$2"/"$z"'_R2_001.fq.gz' trimout/'trmPr_'"$z"'_R1_001.fq.gz' trimout/'trmS_'"$z"'_R1_001.fq.gz' trimout/'trmPr_'"$z"'_R2_001.fq.gz' trimout/'trmS_'"$z"'_R2_001.fq.gz' ILLUMINACLIP:$4:2:30:10 $cuts
	     mv trimfile.txt trimmersettings.txt
	     exit
        elif [[ "$x" =~ .*_R1_001\.fq$ ]]; 
             then
             z=$(basename $x _R1_001.fq)
             java -jar -Xmx1024m /root/Trimmomatic-0.36/trimmomatic-0.36.jar PE "$2"/"$z"'_R1_001.fq' "$2"/"$z"'_R2_001.fq' trimout/'trmPr_'"$z"'_R1_001.fq' trimout/'trmS_'"$z"'_R1_001.fq' trimout/'trmPr_'"$z"'_R2_001.fq' trimout/'trmS_'"$z"'_R2_001.fq' ILLUMINACLIP:$4:2:30:10 $cuts
	     pigz -p 4 trimout/'trmPr_'"$z"'_R1_001.fq' 
	     pigz -p 4 trimout/'trmS_'"$z"'_R1_001.fq' 
	     pigz -p 4 trimout/'trmPr_'"$z"'_R2_001.fq' 
	     pigz -p 4 trimout/'trmS_'"$z"'_R2_001.fq'
             mv trimfile.txt trimmersettings.txt
             exit        
        elif [[ "$x" =~ .*_R1_001\.fastq.gz$ ]]; 
             then
             z=$(basename $x _R1_001.fastq.gz)
             java -jar -Xmx1024m /root/Trimmomatic-0.36/trimmomatic-0.36.jar PE "$2"/"$z"'_R1_001.fastq.gz' "$2"/"$z"'_R2_001.fastq.gz' trimout/'trmPr_'"$z"'_R1_001.fastq.gz' trimout/'trmS_'"$z"'_R1_001.fastq.gz' trimout/'trmPr_'"$z"'_R2_001.fastq.gz' trimout/'trmS_'"$z"'_R2_001.fastq.gz' ILLUMINACLIP:$4:2:30:10 $cuts
	     mv trimfile.txt trimmersettings.txt
             exit
        elif [[ "$x" =~ .*_R1_001\.fastq$ ]]; 
             then
             z=$(basename $x _R1_001.fastq)
             java -jar -Xmx1024m /root/Trimmomatic-0.36/trimmomatic-0.36.jar PE "$2"/"$z"'_R1_001.fastq' "$2"/"$z"'_R2_001.fastq' trimout/'trmPr_'"$z"'_R1_001.fastq' trimout/'trmS_'"$z"'_R1_001.fastq' trimout/'trmPr_'"$z"'_R2_001.fastq' trimout/'trmS_'"$z"'_R2_001.fastq' ILLUMINACLIP:$4:2:30:10 $cuts
	     pigz -p 4 trimout/'trmPr_'"$z"'_R1_001.fastq'
             pigz -p 4 trimout/'trmS_'"$z"'_R1_001.fastq'
             pigz -p 4 trimout/'trmPr_'"$z"'_R2_001.fastq'
             pigz -p 4 trimout/'trmS_'"$z"'_R2_001.fastq'
             mv trimfile.txt trimmersettings.txt
             exit
	else
	     rm -r trimout trimfile.txt
             echo "The input folder - $2 should contain reads that have extension _R1_001.fq.gz/_R2_001.fq.gz or _R1_001.fastq.gz/_R1_001.fastq.gz or _R1_001.fq/_R2_001.fq or _R1_001.fastq/_R2_001.fastq
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
		if [[ "$y" =~ \.gz$ ]]; then
                	java -jar -Xmx1024m /root/Trimmomatic-0.36/trimmomatic-0.36.jar SE $x trimout/$z ILLUMINACLIP:$4:2:30:10 $cuts
                else
                        java -jar -Xmx1024m /root/Trimmomatic-0.36/trimmomatic-0.36.jar SE $x trimout/$z ILLUMINACLIP:$4:2:30:10 $cuts
                        pigz-2.3.4/pigz -p 4 trimout/$z
		fi
            done
        mv trimfile.txt trimmersettings.txt
fi
