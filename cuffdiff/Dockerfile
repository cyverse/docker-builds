# Cuffdiff 2.2.1
#
FROM gcc:5.1 
MAINTAINER Roger Barthelson LABEL Description="This image is used for running Cuffdiff"
RUN git clone https://github.com/samtools/htslib.git && cd htslib && make && make install 
RUN git clone https://github.com/samtools/samtools.git && cd samtools && make && make install
RUN wget http://cole-trapnell-lab.github.io/cufflinks/assets/downloads/cufflinks-2.2.1.Linux_x86_64.tar.gz && tar xvzf cufflinks-2.2.1.Linux_x86_64.tar.gz && cd cufflinks-2.2.1.Linux_x86_64 && cp cuf* /usr/bin && cp gffread /usr/bin && cp gtf_to_sam /usr/bin
RUN wget http://de.iplantcollaborative.org/dl/d/DA72F9C3-8005-4F9C-AFB8-D5165427BC87/munge_ids.pl && wget http://de.iplantcollaborative.org/dl/d/80660887-0DC8-4D25-8840-73D1BED6B90A/cuffdiff_16.pl
RUN wget http://de.iplantcollaborative.org/dl/d/07B1C95B-839C-472D-9DD1-08F840462E95/cuffdiff_sort.sh && wget http://de.iplantcollaborative.org/dl/d/83A9976C-8DDB-4B7D-B540-52AA0DB3774B/cuffdiff_fix_anno.pl && chmod a+x *.pl && chmod a+x cuffdiff_sort.sh
ENTRYPOINT [ "/cuffdiff_16.pl" ]
