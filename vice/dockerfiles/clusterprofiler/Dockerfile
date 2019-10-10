FROM discoenv/rstudio-nginx:3.5.0
MAINTAINER Sateesh Peri <perisateesh@nevada.unr.edu>
LABEL Description "This Dockerfile is for building clusterprofiler"

RUN apt-get update && apt-get install -y zlib1g-dev libudunits2-dev
RUN Rscript -e 'install.packages("littler", dependencies=TRUE)'
RUN Rscript -e 'source("http://bioconductor.org/biocLite.R"); biocLite("AnnotationDbi"); biocLite("org.Mm.eg.db"); biocLite("clusterProfiler"); biocLite("org.Hs.eg.db"); biocLite("org.Mm.eg.db"); biocLite("org.Mmu.eg.db"); biocLite("org.Pt.eg.db"); biocLite("org.Rn.eg.db"); biocLite("org.Sc.sgd.db"); biocLite("org.Ss.eg.db"); biocLite("org.Xl.eg.db"); biocLite("org.Ag.eg.db"); biocLite("org.At.tair.db"); biocLite("org.Bt.eg.db"); biocLite("org.Ce.eg.db"); biocLite("org.Cf.eg.db"); biocLite("org.Dm.eg.db"); biocLite("org.Dr.eg.db"); biocLite("org.EcK12.eg.db"); biocLite("org.EcSakai.eg.db"); biocLite("org.Gg.eg.db");'
RUN wget -qO - https://packages.irods.org/irods-signing-key.asc | apt-key add - \
    && echo "deb [arch=amd64] https://packages.irods.org/apt/ xenial main" > /etc/apt/sources.list.d/renci-irods.list \
    && apt-get update \
    && apt-get install -y irods-icommands \
    && apt-get clean \
    && rm -rf /usr/lib/apt/lists/* 
