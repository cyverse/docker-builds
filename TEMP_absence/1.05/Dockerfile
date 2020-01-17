FROM ubuntu:18.04
MAINTAINER Reetu Tuteja <reetututeja@cyverse.org>
LABEL Description="This image is used for running TEMP"

RUN apt-get update && apt-get install --yes build-essential gcc-multilib apt-utils zlib1g-dev wget samtools bedtools bioperl git
RUN wget https://sourceforge.net/projects/bio-bwa/files/bwa-0.7.17.tar.bz2/download
RUN tar xvf download
WORKDIR /bwa-0.7.17
RUN make
ENV PATH="/bwa-0.7.17:${PATH}"
WORKDIR /

RUN wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/twoBitToFa
RUN chmod +x ./twoBitToFa
RUN mv twoBitToFa /usr/local/bin

RUN git clone https://github.com/JialiUMassWengLab/TEMP.git
WORKDIR /TEMP/scripts
RUN chmod +x ./TEMP_Absence.sh
ENV PATH="/TEMP/scripts:${PATH}"
WORKDIR /

ENTRYPOINT ["TEMP_Absence.sh"]

