FROM java:7

MAINTAINER Roger Barthelson <rogerab@cyverse.org>



RUN apt-get update

RUN apt-get upgrade -y



RUN wget -P /root http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/Trimmomatic-0.33.zip

RUN unzip /root/Trimmomatic-0.33.zip -d /root

RUN wget https://github.com/azroger/iPlant-basic-scripts/archive/v0.0.1.tar.gz && tar xvzf v0.0.1.tar.gz && chmod a+x iPlant-basic-scripts-0.0.1/trimscript.sh && cp iPlant-basic-scripts-0.0.1/trimscript.sh /usr/bin

ENV PATH "/root/Trimmomatic-0.33/adapters:$PATH"



ENTRYPOINT ["trimscript.sh"]