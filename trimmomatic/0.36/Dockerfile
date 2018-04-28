FROM java:7
MAINTAINER Upendra Devisetty <upendra@cyverse.org>
LABEL Description "This Dockerfile is for trimmomatic-0.36 version"

RUN apt-get update
RUN apt-get upgrade -y
# Trimmomatic
RUN wget -P /root http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/Trimmomatic-0.36.zip
RUN unzip /root/Trimmomatic-0.36.zip -d /root
ADD trimscript.wrapper.sh /usr/bin
RUN chmod a+x /usr/bin/trimscript.wrapper.sh
ENV PATH "/root/Trimmomatic-0.36/adapters:$PATH"
# pigz
RUN wget http://www.zlib.net/pigz/pigz-2.3.4.tar.gz
RUN tar xvf pigz-2.3.4.tar.gz
WORKDIR pigz-2.3.4
RUN apt-get install -y make g++ zlib1g-dev
RUN make
RUN cp pigz /usr/bin

ENTRYPOINT ["trimscript.wrapper.sh"]
