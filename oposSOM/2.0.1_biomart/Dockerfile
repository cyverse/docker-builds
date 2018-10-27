FROM rocker/rstudio:3.4.4
MAINTAINER Upendra Devisetty <upendra@cyverse.org>
LABEL version="2.0.1" description="This image is for opossom-2.0.1"

# Run updates 
RUN apt-get update && apt-get install -y libxml2-dev

# Install packages
RUN Rscript -e 'install.packages(c("XML", "igraph", "fastICA", "tsne", "scatterplot3d", "pixmap", "fdrtool", "ape", "RcppParallel", "dplyr"));'
RUN Rscript -e 'source("http://bioconductor.org/biocLite.R"); biocLite("stringi");'
RUN Rscript -e 'install.packages(c("RCurl", "httr"));'
RUN Rscript -e 'source("http://bioconductor.org/biocLite.R"); biocLite("biomaRt");'
RUN Rscript -e 'install.packages("getopt");'

# Copy files to the container
COPY opossom_wrapper.R oposSOM-plants-ensembl_gene_id.tar.gz /

# Install opossom
RUN R CMD INSTALL --build oposSOM-plants-ensembl_gene_id.tar.gz

# change permissions to the wrapper script
ENV BINPATH /usr/bin
RUN chmod +x /opossom_wrapper.R && cp /opossom_wrapper.R $BINPATH

# Add biomart file inside the docker image
COPY biomart_files.txt /

ENTRYPOINT ["opossom_wrapper.R"]
