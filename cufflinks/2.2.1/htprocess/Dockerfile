# Cufflinks 2.21
#
FROM ubuntu:14.04 
MAINTAINER Roger Barthelson <rogerab@email.arizona.edu>
LABEL Description="This image is used for running the HTP_cufflinks-2.2.1 image"
#RUN git clone https://github.com/samtools/htslib.git && cd htslib && make && make install 
#RUN git clone https://github.com/samtools/samtools.git && cd samtools && make && make install
RUN apt-get -y update && apt-get -y install samtools && apt-get -y install wget && apt-get -y install python
RUN wget http://cole-trapnell-lab.github.io/cufflinks/assets/downloads/cufflinks-2.2.1.Linux_x86_64.tar.gz && tar xvzf cufflinks-2.2.1.Linux_x86_64.tar.gz && cd cufflinks-2.2.1.Linux_x86_64 && cp cuf* /usr/bin && cp gffread /usr/bin && cp gtf_to_sam /usr/bin
RUN wget https://github.com/azroger/HTProcess-Pipeline/archive/0.0.7.tar.gz && tar xvzf 0.0.7.tar.gz && chmod a+x HTProcess-Pipeline-0.0.7/* && cp HTProcess-Pipeline-0.0.7/* /usr/bin
ENTRYPOINT [ "HTProcess_cufflinks.sh" ] 
