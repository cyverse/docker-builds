FROM ubuntu:18.04
  
RUN ln -snf /usr/share/zoneinfo/Etc/UTC /etc/localtime \
    && echo "Etc/UTC" > /etc/timezone \
    && apt-get update \
    && apt-get upgrade -y \
    && apt-get install python3-pip texlive-latex-extra wget git dvipng -y \
    && rm -rf /var/lib/apt/lists/*
RUN pip3 install biopython==1.70 numpy deap networkx matplotlib jcvi

# liftover
RUN echo 'export PATH=/opt/conda/bin:$PATH' > /etc/profile.d/conda.sh && \
    wget --quiet https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O ~/miniconda.sh && \
    /bin/bash ~/miniconda.sh -b -p /opt/conda && \
    rm ~/miniconda.sh
ENV PATH /opt/conda/bin:$PATH
RUN conda install -c bioconda ucsc-liftover

# ALLPATHS
RUN mkdir code
WORKDIR code
RUN git clone git://github.com/tanghaibao/jcvi.git
ENV PYTHONPATH /code:$PYTHONPATH
