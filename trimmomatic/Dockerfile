FROM java:7
MAINTAINER Andre Mercer <amercer@iplantcollaborative.org>

RUN apt-get update
RUN apt-get upgrade -y

RUN wget -P /root http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/Trimmomatic-0.33.zip
RUN unzip /root/Trimmomatic-0.33.zip -d /root

ENV PATH "/root/Trimmomatic-0.33/adapters:$PATH"

ENTRYPOINT ["java", "-Xmx1024m", "-jar", "/root/Trimmomatic-0.33/trimmomatic-0.33.jar"]
