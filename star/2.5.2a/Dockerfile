# STAR 
# VERSION              2.5.2.a
#

FROM      ubuntu:14.04.3 
MAINTAINER Kapeel Chougule <kchougul@cshl.edu>

LABEL Description="This image is used for running STAR RNA seq alignemnt tool"
RUN apt-get update && apt-get install -y build-essential wget git curl unzip zlib1g-dev
ADD STAR_align.pl /usr/bin/ 
RUN [ "chmod", "+x", "/usr/bin/STAR_align.pl" ]
RUN wget https://github.com/alexdobin/STAR/archive/2.5.2a.tar.gz
RUN tar zxvf 2.5.2a.tar.gz

WORKDIR /STAR-2.5.2a
RUN make STAR
RUN cp source/STAR /usr/bin/ 

WORKDIR /

ENTRYPOINT ["STAR_align.pl"]

# Build
# sudo docker build -t"=ubuntu/star:2.5.2.a" .
# Test run
# sudo docker run --rm -v $(pwd):/working-dir -w /working-dir ubuntu/star:2.5.2.a -file_query reads_1.fastq -user_database sorghum_bicolor.genome.fasta -user_annotation Sorghum_bicolor.Sorbi1.31.chr.gtf -file_type "SE"
# sudo docker run --rm -v $(pwd):/working-dir -w /working-dir ubuntu/star:2.5.2.a -file_query test.fq -user_database test.fasta -user_annotation test.gtf -file_type "SE"
