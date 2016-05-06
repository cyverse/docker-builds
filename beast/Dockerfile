# Pull base image
FROM java:7

MAINTAINER Andre Mercer <amercer@iplantcollaborative.org>

# Run updates
RUN apt-get update
RUN apt-get upgrade -y
RUN apt-get install wget -y

# Pull down BEAST v1.8.3 
RUN wget -P /root https://github.com/beast-dev/beast-mcmc/releases/download/v1.8.3/BEASTv1.8.3.tgz
RUN tar zxf /root/BEASTv1.8.3.tgz -C /root

# Entrypoint for beast.jar. Will be the default entrypoint for container
ENTRYPOINT ["java", "-Xmx1024m", "-jar", "/root/BEASTv1.8.3/lib/beast.jar"]

# Entrypoint for beauti.jar. Will be added through the Belphegor
#ENTRYPOINT ["java", "-Xmx1024m", "-jar", "/root/BEASTv1.8.3/lib/beauti.jar"]
