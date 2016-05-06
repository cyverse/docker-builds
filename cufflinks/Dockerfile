## Dockerfile
FROM ubuntu:14.04.3
MAINTAINER Upendra Devisetty
LABEL Description="cufflinks-2.2.1.Linux update"

# Install all the updates and download dependencies
RUN apt-get update && apt-get install -y git wget

# Download the cufflinks-2.2.1
RUN wget -O- http://cole-trapnell-lab.github.io/cufflinks/assets/downloads/cufflinks-2.2.1.Linux_x86_64.tar.gz | tar xzvf -

# Clone the cufflinks-2.2.1 wrapper script
RUN git clone https://github.com/upendrak/cufflinks-2.2.1.git 

# Change the permissions and the path for the wrapper script
RUN chmod +x /cufflinks-2.2.1/cufflinks-2.2.1.pl && cp /cufflinks-2.2.1/cufflinks-2.2.1.pl /usr/bin 
ENV PATH /cufflinks-2.2.1.Linux_x86_64/:$PATH

# Entrypoint
ENTRYPOINT ["/usr/bin/cufflinks-2.2.1.pl"]
