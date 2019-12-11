## Dockerfile
FROM ubuntu:16.04
MAINTAINER Amanda Cooksey	
LABEL Description="parse BLAST Bpo--part of the Ortho MCL pipeline"

# Install all the updates and download dependencies
RUN apt-get update
RUN apt-get install --yes \
 apt-utils \
 build-essential \
 libexpat1-dev \
 perl 

RUN cpan App::cpanminus

RUN cpanm --force Bio::Perl XML::Parser XML::DOM Bio::SearchIO

#RUN cpanm -n \
#       Dist::Zilla \
#       Dist::Zilla::Plugin::Git \
#       Dist::Zilla::PluginBundle::BioPerl
 
ADD parseBlastBpo.pl /usr/bin

# Change the permissions and the path for the script
RUN chmod +x /usr/bin/parseBlastBpo.pl

# Entrypoint
ENTRYPOINT ["/usr/bin/parseBlastBpo.pl"]

RUN mkdir /work-dir

WORKDIR /work-dir

