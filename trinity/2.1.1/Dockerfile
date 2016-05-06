# Pull base image
FROM python:2.7

MAINTAINER Andre Mercer <amercer@iplantcollaborative.org>
# Run updates
RUN apt-get update
RUN apt-get upgrade -y
RUN apt-get install wget -y

# Install dependencies
RUN apt-get install build-essential -y
RUN apt-get install samtools -y
RUN apt-get install unzip

RUN wget -P /root 'http://iweb.dl.sourceforge.net/project/bowtie-bio/bowtie/1.1.2/bowtie-1.1.2-linux-x86_64.zip'
RUN unzip /root/bowtie-1.1.2-linux-x86_64.zip -d /root
ENV PATH /root/bowtie-1.1.2:$PATH

RUN wget https://github.com/trinityrnaseq/trinityrnaseq/archive/v2.1.1.zip
RUN unzip v2.1.1.zip -d /root
RUN make -C /root/trinityrnaseq-2.1.1
RUN make plugins -C /root/trinityrnaseq-2.1.1

ENTRYPOINT ["/root/trinityrnaseq-2.1.1/Trinity"]
