FROM ubuntu:24.04
########################################################################################
RUN apt-get update
# Python:
RUN apt-get install -y python3 python3-pip python3-venv
ENV VIRTUAL_ENV=/opt/venv
RUN python3 -m venv $VIRTUAL_ENV
ENV PATH="$VIRTUAL_ENV/bin:$PATH"
# For GeneBlocks:
RUN apt-get install -y ncbi-blast+
# weighted-levenshtein:
RUN apt-get install -y gcc
RUN pip install --no-cache-dir weighted-levenshtein
########################################################################################
# http://www.htslib.org/download/
RUN apt-get install -y libncurses-dev libbz2-dev liblzma-dev
ADD https://github.com/samtools/samtools/releases/download/1.21/samtools-1.21.tar.bz2 .
RUN tar xjvf samtools-1.21.tar.bz2
RUN cd samtools-1.21 \
    && ./configure --prefix=/opt/samtools \
    && make \
    && make install
ENV PATH="/opt/samtools/bin:$PATH"
########################################################################################
ADD https://github.com/samtools/bcftools/releases/download/1.21/bcftools-1.21.tar.bz2 .
RUN tar xjvf bcftools-1.21.tar.bz2
RUN cd bcftools-1.21 \
    && ./configure --prefix=/opt/bcftools \
    && make \
    && make install
ENV PATH="/opt/bcftools/bin:$PATH"
########################################################################################
ADD https://github.com/samtools/htslib/releases/download/1.21/htslib-1.21.tar.bz2 .
RUN tar xjvf htslib-1.21.tar.bz2
RUN cd htslib-1.21 \
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
ADD https://github.com/lh3/minimap2/releases/download/v2.28/minimap2-2.28_x64-linux.tar.bz2 .
RUN tar -jxvf minimap2-2.28_x64-linux.tar.bz2 \
    && mv minimap2-2.28_x64-linux/ /opt/minimap2-2.28_x64-linux
ENV PATH="/opt/minimap2-2.28_x64-linux:$PATH"
########################################################################################
# Java:
RUN apt-get install -y default-jre
ENV JAVA_HOME=/usr/lib/jvm/default-java/
########################################################################################
# Canu:
# Install from source...
# RUN apt-get install -y libcurl4-gnutls-dev libpcap-dev libssl-dev
# RUN apt-get install -y git-all
# RUN git clone https://github.com/marbl/canu.git \
#     && cd canu/src \
#     && make -j
# ENV PATH="/canu/build/bin:$PATH"
########
# Simple canu install from release (once the patch release is made https://github.com/marbl/canu/issues/2359#issuecomment-2631317846)
ADD https://github.com/marbl/canu/releases/download/v2.2/canu-2.2.Linux-amd64.tar.xz .
RUN tar -xJf canu-2.2.Linux-amd64.tar.xz \
    && mv canu-2.2/ /opt/canu-2.2
ENV PATH="/opt/canu-2.2/bin:$PATH"
# Latest release:
# ADD https://github.com/marbl/canu/releases/download/v2.3/canu-2.3.Linux-amd64.tar.xz .
# RUN tar -xJf canu-2.3.Linux-amd64.tar.xz \
#     && mv canu-2.3/ /opt/canu-2.3
# ENV PATH="/opt/canu-2.3/bin:$PATH"
########################################################################################
# For DNA Features Viewer:
RUN pip install --no-cache-dir packaging
# Python packages:
RUN pip install --no-cache-dir networkx==3.4.2 python-Levenshtein==0.26.1 pypugjs==5.11.0 jinja2==3.1.4 weasyprint==64.0 beautifulsoup4==4.12.3 Markdown==3.7
RUN pip install --no-cache-dir matplotlib==3.10.0 numpy==2.2.0 pandas==2.2.3 weighted_levenshtein==0.2.2 portion==2.6.0
RUN pip install --no-cache-dir biopython==1.84
RUN pip install --no-cache-dir dna_features_viewer==3.1.3 geneblocks==1.2.3 cyvcf2==0.31.1 pdf_reports==0.3.7
RUN pip install --no-cache-dir NanoPlot==1.44.1 nanofilt==2.8.0
RUN apt-get install -y git-all
RUN pip install ediacara@git+https://github.com/Edinburgh-Genome-Foundry/Ediacara.git
RUN pip install --no-cache-dir pypdf==5.3.0
RUN pip install --no-cache-dir pydyf==0.11.0
########################################################################################
# Opencontainer labels
LABEL org.opencontainers.image.title="sequeduct"
LABEL org.opencontainers.image.description="Sequencing analysis pipeline for validating DNA assembly constructs"
LABEL org.opencontainers.image.url="https://github.com/edinburgh-genome-foundry/sequeduct"
LABEL org.opencontainers.image.documentation="https://github.com/edinburgh-genome-foundry/sequeduct"
LABEL org.opencontainers.image.source="https://github.com/edinburgh-genome-foundry/sequeduct/"
LABEL org.opencontainers.image.vendor="edinburgh-genome-foundry"
LABEL org.opencontainers.image.authors="Peter Vegh"
LABEL org.opencontainers.image.revision="v0.4.4"
