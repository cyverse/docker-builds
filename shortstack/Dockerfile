FROM centos:7

MAINTAINER Mike Axtell mja18@psu.edu

#install make
RUN yum -y install make

#install gcc and gcc-c++
RUN yum -y install gcc gcc-c++

#install unzip and bzip2
RUN yum -y install unzip bzip2

# install zlib
ADD http://zlib.net/zlib-1.2.8.tar.gz /
RUN tar xvf zlib-1.2.8.tar.gz && cd zlib-1.2.8 && ./configure && make && \
make install

#install curses
ADD http://ftp.gnu.org/pub/gnu/ncurses/ncurses-5.4.tar.gz /
RUN tar xvf ncurses-5.4.tar.gz && cd ncurses-5.4 && ./configure && make && \
make install

# download and install samtools
ADD https://github.com/samtools/samtools/releases/download/1.2/samtools-1.2.tar.bz2 /
RUN tar xvf samtools-1.2.tar.bz2 && cd samtools-1.2 && make && make install

#install perl
RUN yum -y install perl-5.16.3

#install csh
RUN yum -y install csh

#install gd
RUN yum -y install gd

#install RNAfold
ADD http://www.tbi.univie.ac.at/RNA/download.php?id=viennarna-2.1.9-x86_64-rpm viennarna-2.1.9.rpm
RUN rpm -Uvh viennarna-2.1.9.rpm

#install bowtie
ADD http://downloads.sourceforge.net/project/bowtie-bio/bowtie/1.1.1/bowtie-1.1.1-linux-x86_64.zip?r=http%3A%2F%2Fsourceforge.net%2Fprojects%2Fbowtie-bio%2Ffiles%2Fbowtie%2F1.1.1%2F&ts=1440521017&use_mirror=skylineservers bowtie-1.1.1.zip
RUN unzip bowtie-1.1.1.zip
ENV PATH /bowtie-1.1.1:$PATH

#install latest ShortStack from github
ADD https://raw.githubusercontent.com/MikeAxtell/ShortStack/master/ShortStack /usr/bin/ShortStack

RUN chmod 755 /usr/bin/ShortStack
