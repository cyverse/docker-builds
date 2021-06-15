## Dockerfile
FROM ubuntu:18.04
MAINTAINER Amanda Cooksey	
LABEL Description="fastq-multx"

# Install all the updates and download dependencies
RUN apt-get update && \
    apt-get install -y \
    git \
    wget \
    bzip2

RUN echo 'export PATH=/opt/conda/bin:$PATH' > /etc/profile.d/conda.sh && \
    wget --quiet https://repo.continuum.io/miniconda/Miniconda2-4.0.5-Linux-x86_64.sh -O ~/miniconda.sh && \
    /bin/bash ~/miniconda.sh -b -p /opt/conda && \
    rm ~/miniconda.sh


# give write permissions to conda folder
RUN chmod 777 -R /opt/conda/

ENV PATH=$PATH:/opt/conda/bin

RUN conda config --add channels bioconda

RUN conda upgrade conda

# add blast 2.7.1--i added this
RUN conda install -c conda-forge -c bioconda fastq-multx==1.4.0

ENV PATH /masterscript.sh/:$PATH

# add masterscript.sh
ADD masterscript.sh /usr/bin

# Change the permissions and the path for the wrapper script
RUN chmod +x /usr/bin/masterscript.sh

RUN mkdir /work-dir

# Entrypoint
ENTRYPOINT ["/usr/bin/masterscript.sh"]


# Add path to working directory
WORKDIR /work-dir
