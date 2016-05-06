# Fastqc-plus
#
FROM ubuntu 
MAINTAINER Roger Barthelson LABEL Description="This image is used for running FastQC and creating an output html file"
RUN apt-get -y update
RUN apt-get -y install wget && apt-get -y install unzip && apt-get -y install make
RUN wget http://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v0.10.1.zip && unzip fastqc_v0.10.1.zip && chmod a+x FastQC/fastqc && cp -r FastQC/* /usr/bin
RUN apt-get -y install tzdata-java openjdk-7-jre-headless openjdk-7-jre default-jre
RUN wget -O - http://cpanmin.us | perl - App::cpanminus && cpanm File::Slurp
RUN wget https://github.com/azroger/iPlant-basic-scripts/archive/master.zip && unzip master.zip
RUN chmod a+x iPlant-basic-scripts-master/embed_images.pl && chmod a+x iPlant-basic-scripts-master/fastqc_plus.sh
RUN cp iPlant-basic-scripts-master/fastqc_plus.sh /usr/bin && cp iPlant-basic-scripts-master/embed_images.pl /usr/bin
ENTRYPOINT [ "fastqc_plus.sh" ]
