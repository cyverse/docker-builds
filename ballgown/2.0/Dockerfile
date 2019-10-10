FROM kapeel/ballgown-r-package
MAINTAINER Upendra Devisetty <upendra@cyverse.org>
LABEL Description="This image is used for Ballgown RNAseq differential expression analysis"
# Install dependencies
RUN apt-get update 
RUN Rscript -e 'install.packages(c("dplyr","ggplot2", "reshape2", "cowplot"), dependencies = TRUE);'
ADD ballgown.R /usr/bin/
RUN [ "chmod", "+x",  "/usr/bin/ballgown.R" ]
ENTRYPOINT ["ballgown.R"]
