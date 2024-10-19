FROM python:3.9
########################################################################################
RUN apt-get update
# For GeneBlocks:
RUN apt-get install -y ncbi-blast+
# weighted-levenshtein:
RUN apt-get install -y gcc
RUN pip install --no-cache-dir weighted-levenshtein
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
# Canu:
ADD https://github.com/marbl/canu/releases/download/v2.2/canu-2.2.Linux-amd64.tar.xz .
RUN tar -xJf canu-2.2.Linux-amd64.tar.xz \
    && mv canu-2.2/ /opt/canu-2.2
ENV PATH="/opt/canu-2.2/bin:$PATH"
########################################################################################
# Java:
RUN apt-get install -y default-jre
ENV JAVA_HOME /usr/lib/jvm/default-java/
########################################################################################
# For DNA Features Viewer:
RUN pip install --no-cache-dir packaging
# Python packages:
RUN pip install --no-cache-dir networkx==3.1 python-Levenshtein==0.21.1 pypugjs==5.9.12 jinja2==3.1.2 weasyprint==59.0 beautifulsoup4==4.12.2 Markdown==3.4.4
RUN pip install --no-cache-dir matplotlib==3.7.2 numpy==1.25.2 pandas==2.0.3 weighted_levenshtein==0.2.2 portion==2.4.1
RUN pip install --no-cache-dir biopython==1.78
RUN pip install --no-cache-dir dna_features_viewer==3.1.2 geneblocks==1.2.3 cyvcf2==0.30.22 pdf_reports==0.3.5
RUN pip install --no-cache-dir geneblocks==1.2.3
RUN pip install --no-cache-dir NanoPlot==1.38.1 nanofilt==2.8.0
RUN pip install ediacara@git+https://github.com/Edinburgh-Genome-Foundry/Ediacara.git
RUN pip install --no-cache-dir pypdf==4.2.0
RUN pip install --no-cache-dir pydyf==0.10.0
########################################################################################
# Opencontainer labels
LABEL org.opencontainers.image.title="sequeduct"
LABEL org.opencontainers.image.description="Sequencing analysis pipeline for validating DNA assembly constructs"
LABEL org.opencontainers.image.url="https://github.com/edinburgh-genome-foundry/sequeduct"
LABEL org.opencontainers.image.documentation="https://github.com/edinburgh-genome-foundry/sequeduct"
LABEL org.opencontainers.image.source="https://github.com/edinburgh-genome-foundry/sequeduct/"
LABEL org.opencontainers.image.vendor="edinburgh-genome-foundry"
LABEL org.opencontainers.image.authors="Peter Vegh"
LABEL org.opencontainers.image.revision="v0.3.2"
