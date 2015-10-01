FROM perl:latest

MAINTAINER Ken Youens-Clark <kyclark@email.arizona.edu>

RUN apt-get update && apt-get install libdb-dev -y

RUN cpanm --force Capture::Tiny

RUN cpanm --force BioPerl

RUN cpanm File::Which

COPY wrapper_phage_contigs_sorter_iPlant.pl /usr/local/bin/

COPY Scripts /usr/local/bin/Scripts/

ADD https://github.com/iPlantCollaborativeOpenSource/docker-builds/releases/download/virsorter/virsorter.bin.tgz ./

RUN tar -xz -C /usr/local -f virsorter.bin.tgz && rm virsorter.bin.tgz

ENTRYPOINT ["wrapper_phage_contigs_sorter_iPlant.pl"]

CMD ["-h"]
