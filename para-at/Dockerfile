#modeled on https://github.com/bgruening/docker-recipes/blob/master/coprarna/Dockerfile
FROM ubuntu:14.04
MAINTAINER Jeremy DeBarry   jdebarry@cyverse.org
#Get paraat.pl  and epal2nal.pl code and add into $PATH to make it executable anywhere, then make it executable
ADD https://raw.githubusercontent.com/jdebarry/paraat/master/ParaAT2.0/ParaAT.pl /usr/bin/
ADD https://raw.githubusercontent.com/jdebarry/paraat/master/ParaAT2.0/Epal2nal.pl /usr/bin/
RUN [ "chmod", "+x",  "/usr/bin/ParaAT.pl" ]
RUN [ "chmod", "+x",  "/usr/bin/Epal2nal.pl" ]
#Installing aligners, renaming clustalw executable so paraat.pl will use it - this is a bandaid but it works
RUN echo "deb http://archive.ubuntu.com/ubuntu trusty multiverse" >> /etc/apt/sources.list
RUN DEBIAN_FRONTEND=noninteractive apt-get -qq update
RUN DEBIAN_FRONTEND=noninteractive apt-get install --no-install-recommends -y clustalw
RUN DEBIAN_FRONTEND=noninteractive apt-get install --no-install-recommends -y mafft
RUN DEBIAN_FRONTEND=noninteractive apt-get install --no-install-recommends -y muscle
RUN DEBIAN_FRONTEND=noninteractive apt-get install --no-install-recommends -y t-coffee
RUN [ "mv", "/usr/bin/clustalw" , "usr/bin/clustalw2" ]
ENTRYPOINT ["ParaAT.pl"]
#Command Line Usage: docker run --rm -v=/Users/jdebarry/Dropbox/Docker/paraat/Development/input/:/data paraat -h /data/test.homologs -n /data/test.cds -a /data/test.pep -p proc -o /data/outputclustal -v -m clustalw2 -f axt
#-v is verbose, -m is for aligner choice: tcoffee (default) | clustalw2 | mafft | muscle, -f is for output format: axt | fasta (default) | paml | codon | clustal, optional)
