#Pull base image
FROM r-base

MAINTAINER Andre Mercer <amercer@iplantcollaborative.org>

# Run updates
RUN apt-get update
RUN apt-get upgrade -y
RUN apt-get install wget -y
RUN apt-get install r-base-dev -y
RUN apt-get install libxml2 -y
RUN apt-get install libxml2-dev -y

RUN Rscript -e 'source("http://bioconductor.org/biocLite.R"); biocLite("DESeq2");'
RUN Rscript -e 'install.packages("getopt");'
RUN Rscript -e 'install.packages("gplots");'

COPY run_DESeq2.R /usr/local/bin/

RUN chmod +x /usr/local/bin/run_DESeq2.R

ENTRYPOINT ["run_DESeq2.R"]
