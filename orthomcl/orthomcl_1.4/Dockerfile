## Dockerfile
FROM ubuntu:16.04
MAINTAINER Amanda Cooksey	
LABEL Description="OrthoMCL version 1.4 uses parsed BLASTp input to cluster homologs based on a Markov Clustering Algorithm."

RUN cat /proc/version

# Install all the updates and download dependencies
RUN apt-get update && apt-get install -y \
    perl \
    make \
    build-essential \
    wget \
    apt-utils \
    libexpat1-dev \
    gcc-multilib

RUN cpan App::cpanminus

RUN cpanm --force Bio::SearchIO 

#RUN cpanm -n \
#       Dist::Zilla \
#       Dist::Zilla::Plugin::Git \
#       Dist::Zilla::PluginBundle::BioPerl


RUN mkdir /work-dir

WORKDIR /work-dir

RUN wget https://orthomcl.org/common/downloads/software/unsupported/v1.4/mcl-02-063.tar.gz && tar -xvzf mcl-02-063.tar.gz #&& \
	wget https://orthomcl.org/common/downloads/software/unsupported/v1.4/ORTHOMCL_V1.4_mcl-02-063.tar && tar -xf ORTHOMCL_V1.4_mcl-02-063.tar

WORKDIR /work-dir/mcl-02-063

RUN  /work-dir/mcl-02-063/configure && \
	make && \
	make install #&& \
	make clean


ADD orthomcl.sh orthomcl.pl /usr/bin/ 
ADD orthomcl_io.pm orthomcl_core.pm /etc/perl/

# Change the permissions and the path for the script
RUN chmod a+x  /usr/bin/orthomcl.pl /usr/bin/orthomcl.sh /usr/local/bin/mcl

WORKDIR /work-dir

# Entrypoint
#ENTRYPOINT ["/usr/bin/orthomcl.sh"]
