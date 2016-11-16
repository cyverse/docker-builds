FROM ubuntu:14.04
MAINTAINER Andre Mercer <amercer@iplantcollaborative.org>
LABEL Description "This Dockerfile is for installing trim-galore:0.4.1"

# Run updates
RUN apt-get clean
RUN apt-get update
RUN apt-get upgrade -y
RUN apt-get install -y wget unzip python openjdk-7-jre

# Install trim-galore
RUN wget http://www.bioinformatics.babraham.ac.uk/projects/trim_galore/trim_galore_v0.4.1.zip && \
        unzip trim_galore_v0.4.1.zip
ENV PATH /trim_galore_zip:$PATH

# Install cutadapt
RUN wget https://github.com/ComputationalSystemsBiology/EoulsanDockerFiles/raw/master/TrimAdapt/cutadapt-1.8.1.tar.gz && \
        tar -xzf cutadapt-1.8.1.tar.gz
ENV PATH /cutadapt-1.8.1/bin:$PATH

# Install fastqc
RUN wget http://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v0.11.5.zip && \
        unzip fastqc_v0.11.5.zip && \
        cd FastQC && \
        chmod a+x fastqc
ENV PATH /FastQC:$PATH

ENTRYPOINT ["trim_galore"]
CMD ["-h"]
