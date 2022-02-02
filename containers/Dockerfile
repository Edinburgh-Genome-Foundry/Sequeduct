FROM python:3.6
########################################################################################
RUN apt-get update
# python-Levenshtein:
RUN apt-get install -y gcc
RUN pip install --no-cache-dir weighted-levenshtein
# GeneBlocks:
RUN pip install --no-cache-dir biopython==1.78
RUN apt-get install -y ncbi-blast+
RUN pip install --no-cache-dir geneblocks
# Python packages
RUN pip install --no-cache-dir NanoPlot==1.38.1 nanofilt==2.8.0
########################################################################################
# http://www.htslib.org/download/
ADD https://github.com/samtools/samtools/releases/download/1.14/samtools-1.14.tar.bz2 .
RUN tar xjvf samtools-1.14.tar.bz2
RUN cd samtools-1.14 \
    && ./configure --prefix=/opt/samtools \
    && make \
    && make install
ENV PATH="/opt/samtools/bin:$PATH"
########################################################################################
ADD https://github.com/samtools/bcftools/releases/download/1.14/bcftools-1.14.tar.bz2 .
RUN tar xjvf bcftools-1.14.tar.bz2
RUN cd bcftools-1.14 \
    && ./configure --prefix=/opt/bcftools \
    && make \
    && make install
ENV PATH="/opt/bcftools/bin:$PATH"
########################################################################################
ADD https://github.com/samtools/htslib/releases/download/1.14/htslib-1.14.tar.bz2 .
RUN tar xjvf htslib-1.14.tar.bz2
RUN cd htslib-1.14 \
    && ./configure --prefix=/opt/htslib \
    && make \
    && make install
ENV PATH="/opt/htslib/bin:$PATH"
########################################################################################
ADD https://github.com/freebayes/freebayes/releases/download/v1.3.6/freebayes-1.3.6-linux-amd64-static.gz .
RUN gunzip freebayes-1.3.6-linux-amd64-static.gz \
    && chmod +x freebayes-1.3.6-linux-amd64-static \
    && mkdir --parents /opt/freebayes/bin \
    && mv freebayes-1.3.6-linux-amd64-static /opt/freebayes/bin/freebayes
ENV PATH="/opt/freebayes/bin:$PATH"
########################################################################################
ADD https://github.com/lh3/minimap2/releases/download/v2.24/minimap2-2.24_x64-linux.tar.bz2 .
RUN tar -jxvf minimap2-2.24_x64-linux.tar.bz2 \
    && mv minimap2-2.24_x64-linux/ /opt/minimap2-2.24_x64-linux
ENV PATH="/opt/minimap2-2.24_x64-linux:$PATH"
########################################################################################
# For DNA Features Viewer:
RUN pip install packaging
########################################################################################
COPY Ediacara /opt/Ediacara
RUN pip install --no-cache-dir /opt/Ediacara
