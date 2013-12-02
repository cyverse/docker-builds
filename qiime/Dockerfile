FROM ubuntu
ENV DEBIAN_FRONTEND noninteractive
# Work around initramfs-tools running on kernel 'upgrade': <http://bugs.debian.org/cgi-bin/bugreport.cgi?bug=594189>
ENV INITRD No

# Work around initscripts trying to mess with /dev/shm: <https://bugs.launchpad.net/launchpad/+bug/974584>
# Used by our `src/ischroot` binary to behave in our custom way, to always say we are in a chroot.
ENV FAKE_CHROOT 1
RUN mv /usr/bin/ischroot /usr/bin/ischroot.original
ADD src/ischroot /usr/bin/ischroot

# Configure no init scripts to run on package updates.
ADD src/policy-rc.d /usr/sbin/policy-rc.d

RUN echo "deb http://us.archive.ubuntu.com/ubuntu/ precise universe" >> /etc/apt/sources.list
RUN echo "deb http://us.archive.ubuntu.com/ubuntu/ precise multiverse" >> /etc/apt/sources.list
RUN apt-get update
RUN apt-get --force-yes -y install python-dev libncurses5-dev libssl-dev libzmq-dev libgsl0-dev openjdk-6-jdk libxml2 libxslt1.1 libxslt1-dev ant git subversion build-essential zlib1g-dev libpng12-dev libfreetype6-dev mpich2 libreadline-dev gfortran unzip libmysqlclient18 libmysqlclient-dev ghc sqlite3 libsqlite3-dev git python
RUN git clone git://github.com/qiime/qiime-deploy.git /usr/local/qiime-deploy
RUN git clone git://github.com/qiime/qiime-deploy-conf.git /usr/local/qiime-deploy-conf
RUN python /usr/local/qiime-deploy/qiime-deploy.py /usr/local/qiime_software/ -f /usr/local/qiime-deploy-conf/qiime-1.7.0-dev/qiime.conf --force-remove-failed-dirs
