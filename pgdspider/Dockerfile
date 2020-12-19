FROM tswetnam/xpra:bionic

MAINTAINER Amanda Cooksey <amandamcooksey@gmail.com>

USER root
RUN export DEBIAN_FRONTEND=noninteractive

WORKDIR /opt

RUN apt-get update && \
    apt-get install -y openjdk-8-jre

RUN wget http://www.cmpg.unibe.ch/software/PGDSpider/PGDSpider_2.1.1.5.zip \
        && unzip PGDSpider_2.1.1.5.zip \
        && rm PGDSpider_2.1.1.5.zip 

ENV PATH=/opt/PGDSpider_2.1.1.5:$PATH
ENV CLASSPATH=/opt/PGDSpider_2.1.1.5/PGDSpider2.jar

WORKDIR /home/user/

USER user

CMD xpra start --bind-tcp=0.0.0.0:9876 --html=on --start-child="java -Xmx1024m -Xms512m -jar /opt/PGDSpider_2.1.1.5/PGDSpider2.jar" --exit-with-children=no --daemon=no --xvfb="/usr/bin/Xvfb +extension Composite -screen 0 1920x1080x24+32 -nolisten tcp -noreset" --pulseaudio=no --notifications=no --bell=no :100  --video-encoders=nvenc
