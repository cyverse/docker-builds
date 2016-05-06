FROM ubuntu:14.04

MAINTAINER Jeremy DeBarry jdebarry@iplantcollaborative.org

#install make
RUN apt-get -qq install make

#Install clustalw modeled on http://www.microhowto.info/howto/perform_an_unattended_installation_of_a_debian_package.html
RUN echo "deb http://archive.ubuntu.com/ubuntu trusty multiverse" >> /etc/apt/sources.list
RUN DEBIAN_FRONTEND=noninteractive apt-get -qq update
RUN DEBIAN_FRONTEND=noninteractive apt-get install --no-install-recommends -y clustalw

#install cpan minus modeled from http://perlmaven.com/install-perl-modules-without-root-rights-on-linux-ubuntu-13-10
RUN apt-get -qq install curl && curl -L http://cpanmin.us | perl - App::cpanminus

# Install BioPerl Run Package modeled from genomicpariscentre/bioperl:latest
RUN cpanm --force CJFIELDS/BioPerl-Run-1.006900.tar.gz

#Install Log4perl modeled on https://hub.docker.com/r/sjackman/assemblers/~/dockerfile/
RUN cpanm Log::Log4perl

#ncbi-blast modeled from https://bcrc.bio.umass.edu/courses/spring2012/micbio/micbio660/content/blast
RUN apt-get -qq install blast2

#Get AAARFv.1.0.1.pl code and add into $PATH to make it executable anywhere, then make it executable
ADD https://github.com/jdebarry/AAARF/raw/master/v1.0.1/AAARFv1.0.1.pl /usr/bin/
RUN [ "chmod", "+x",  "/usr/bin/AAARFv1.0.1.pl" ]

#create directory for running
CMD mkdir data

#Running AAARF - Input file and BLASTdb files must be in same folder and this folder must be used as input folder for command line execution of image
#Command Line Usage of Image: docker run --rm -v=/path/to/host/input/folder:/data image_name --inputFile=/data/InputFileName
ENTRYPOINT ["AAARFv1.0.1.pl"]
