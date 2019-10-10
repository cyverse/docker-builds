FROM gcc

RUN set -e \
      && apt-get -y update \
      && apt-get -y upgrade \
      && apt-get -y install mercurial parallel \
      && apt-get clean

RUN set -e \
      && hg clone http://last.cbrc.jp/last /usr/local/src/last \
      && cd /usr/local/src/last \
      && make \
      && make install

ADD lastal_wrapper_script.sh /usr/bin
RUN chmod +x /usr/bin/lastal_wrapper_script.sh

ENTRYPOINT ["/usr/bin/lastal_wrapper_script.sh"]
