FROM java:7

MAINTAINER Roger Barthelson <rogerab@cyverse.org>

RUN apt-get update

RUN apt-get upgrade -y

RUN wget -P /root http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/Trimmomatic-0.33.zip

RUN unzip /root/Trimmomatic-0.33.zip -d /root

ADD trimscript.wrapper.sh /usr/bin

RUN chmod a+x /usr/bin/trimscript.wrapper.sh

ENV PATH "/root/Trimmomatic-0.33/adapters:$PATH"

ENTRYPOINT ["trimscript.wrapper.sh"]
