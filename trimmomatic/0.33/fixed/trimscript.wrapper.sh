#!/bin/bash

cat $3 | tr '\n' ' ' > trimfile.txt
cuts=`cat trimfile.txt`
mkdir trimout

if [ "$1" == PE ] && [ -z "$4" ]; then
        for x in $2/*_R1_*
            do
                extension=$(echo "$x" | awk -F '.' '{print $(NF-1) "." $(NF)}')

                if [[ "$extension" =~ fastq.gz ]];
                    then
                        z=$(basename $x ".fastq.gz")
                        z2=${z/_R1_/_R2_}
                        java -jar -Xmx1024m /root/Trimmomatic-0.33/trimmomatic-0.33.jar PE "$2"/"$z"'.fastq.gz' "$2"/"$z2"'.fastq.gz' trimout/'trmPr_'"$z"'.fastq.gz' trimout/'trmS_'"$z"'.fastq.gz' trimout/'trmPr_'"$z2"'.fastq.gz' trimout/'trmS_'"$z2"'.fastq.gz' $cuts
                        cat trimfile.txt | tr ' ' '\n' > trimmersettings.txt
                

                elif [[ "$extension" =~ fq.gz ]];
                    then
                        z=$(basename $x ".fq.gz")
                        z2=${z/_R1_/_R2_}
                        java -jar -Xmx1024m /root/Trimmomatic-0.33/trimmomatic-0.33.jar PE "$2"/"$z"'.fq.gz' "$2"/"$z2"'.fq.gz' trimout/'trmPr_'"$z"'.fq.gz' trimout/'trmS_'"$z"'.fq.gz' trimout/'trmPr_'"$z2"'.fq.gz' trimout/'trmS_'"$z2"'.fq.gz' $cuts
                        cat trimfile.txt | tr ' ' '\n' > trimmersettings.txt
                

                elif [[ "$extension" =~ fastq ]];
                    then
                        z=$(basename $x ".fastq")
                        z2=${z/_R1_/_R2_}
                        java -jar -Xmx1024m /root/Trimmomatic-0.33/trimmomatic-0.33.jar PE "$2"/"$z"'.fastq' "$2"/"$z2"'.fastq' trimout/'trmPr_'"$z"'.fastq' trimout/'trmS_'"$z"'.fastq' trimout/'trmPr_'"$z2"'.fastq' trimout/'trmS_'"$z2"'.fastq' $cuts
                        cat trimfile.txt | tr ' ' '\n' > trimmersettings.txt
                

                elif [[ "$extension" =~ fq ]];
                    then
                        z=$(basename $x ".fq")
                        z2=${z/_R1_/_R2_}
                        java -jar -Xmx1024m /root/Trimmomatic-0.33/trimmomatic-0.33.jar PE "$2"/"$z"'.fq' "$2"/"$z2"'.fq' trimout/'trmPr_'"$z"'.fq' trimout/'trmS_'"$z"'.fq' trimout/'trmPr_'"$z2"'.fq' trimout/'trmS_'"$z2"'.fq' $cuts
                        cat trimfile.txt | tr ' ' '\n' > trimmersettings.txt
                    

                fi
            done


elif [ "$1" == PE ] && [ ! -z "$4" ]; then
        for x in $2/*_R1_*
            do
                extension=$(echo "$x" | awk -F '.' '{print $(NF-1) "." $(NF)}')

                if [[ "$extension" =~ fastq.gz ]];
                    then
                        z=$(basename $x ".fastq.gz")
                        z2=${z/_R1_/_R2_}
                        java -jar -Xmx1024m /root/Trimmomatic-0.33/trimmomatic-0.33.jar PE "$2"/"$z"'.fastq.gz' "$2"/"$z2"'.fastq.gz' trimout/'trmPr_'"$z"'.fastq.gz' trimout/'trmS_'"$z"'.fastq.gz' trimout/'trmPr_'"$z2"'.fastq.gz' trimout/'trmS_'"$z2"'.fastq.gz' ILLUMINACLIP:$4:2:30:10 $cuts
                        cat trimfile.txt | tr ' ' '\n' > trimmersettings.txt
                

                elif [[ "$extension" =~ fq.gz ]];
                    then
                        z=$(basename $x ".fq.gz")
                        z2=${z/_R1_/_R2_}
                        java -jar -Xmx1024m /root/Trimmomatic-0.33/trimmomatic-0.33.jar PE "$2"/"$z"'.fq.gz' "$2"/"$z2"'.fq.gz' trimout/'trmPr_'"$z"'.fq.gz' trimout/'trmS_'"$z"'.fq.gz' trimout/'trmPr_'"$z2"'.fq.gz' trimout/'trmS_'"$z2"'.fq.gz' ILLUMINACLIP:$4:2:30:10 $cuts
                        cat trimfile.txt | tr ' ' '\n' > trimmersettings.txt
                

                elif [[ "$extension" =~ fastq ]];
                    then
                        z=$(basename $x ".fastq")
                        z2=${z/_R1_/_R2_}
                        java -jar -Xmx1024m /root/Trimmomatic-0.33/trimmomatic-0.33.jar PE "$2"/"$z"'.fastq' "$2"/"$z2"'.fastq' trimout/'trmPr_'"$z"'.fastq' trimout/'trmS_'"$z"'.fastq' trimout/'trmPr_'"$z2"'.fastq' trimout/'trmS_'"$z2"'.fastq' ILLUMINACLIP:$4:2:30:10 $cuts
                        cat trimfile.txt | tr ' ' '\n' > trimmersettings.txt
                

                elif [[ "$extension" =~ fq ]];
                    then
                        z=$(basename $x ".fq")
                        z2=${z/_R1_/_R2_}
                        java -jar -Xmx1024m /root/Trimmomatic-0.33/trimmomatic-0.33.jar PE "$2"/"$z"'.fq' "$2"/"$z2"'.fq' trimout/'trmPr_'"$z"'.fq' trimout/'trmS_'"$z"'.fq' trimout/'trmPr_'"$z2"'.fq' trimout/'trmS_'"$z2"'.fq' ILLUMINACLIP:$4:2:30:10 $cuts
                        cat trimfile.txt | tr ' ' '\n' > trimmersettings.txt
                    

                fi
            done


elif [ $1 == SE ] && [ -z "$4" ]; then
        for x in $2/*
            do
                y=$(basename $x)
                z='trimS'$y
                java -jar -Xmx1024m /root/Trimmomatic-0.33/trimmomatic-0.33.jar SE $x trimout/$z $cuts
                cat trimfile.txt | tr ' ' '\n' > trimmersettings.txt
        
            done


elif [ $1 == SE ] && [ ! -z "$4" ]; then
        for x in $2/*
            do
                y=$(basename $x)
                z='trimS'$y
                java -jar -Xmx1024m /root/Trimmomatic-0.33/trimmomatic-0.33.jar SE $x trimout/$z ILLUMINACLIP:$4:2:30:10 $cuts
                cat trimfile.txt | tr ' ' '\n' > trimmersettings.txt
        
            done
fi
rm trimfile.txt

