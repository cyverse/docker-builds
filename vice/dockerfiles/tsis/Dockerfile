FROM discoenv/rstudio-nginx:3.5.0
MAINTAINER Upendra Devisetty <upendra@cyverse.org>
LABEL Description "This Dockerfile is for TSIS shiny app"

RUN Rscript -e "install.packages(c('shiny', 'ggplot2', 'zoo', 'shinyFiles', 'gtools'), repos='http://cran.rstudio.com/')"

RUN wget -qO - https://packages.irods.org/irods-signing-key.asc | apt-key add - \
    && echo "deb [arch=amd64] https://packages.irods.org/apt/ xenial main" > /etc/apt/sources.list.d/renci-irods.list \
    && apt-get update \
    && apt-get install -y irods-icommands \
    && apt-get clean \
    && rm -rf /usr/lib/apt/lists/* 
