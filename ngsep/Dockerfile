# Pull base image
FROM java:7

MAINTAINER Andre Mercer <amercer@iplantcollaborative.org>

# Run updates
RUN apt-get update
RUN apt-get upgrade -y

# Copy tool
COPY ngsep_files /root/ngsep

ENTRYPOINT ["java", "-jar", "/root/ngsep/NGSToolsApp_2.1.4.jar"]
