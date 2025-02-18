# For GeneBlocks:
apt-get install -y ncbi-blast+
# weighted-levenshtein:
apt-get install -y gcc
pip install --no-cache-dir weighted-levenshtein
########################################################################################
# http://www.htslib.org/download/
apt-get install -y libncurses-dev libbz2-dev liblzma-dev
# cd to a software dir
wget https://github.com/samtools/samtools/releases/download/1.21/samtools-1.21.tar.bz2 .
tar xjvf samtools-1.21.tar.bz2
cd samtools-1.21 \
    && ./configure \
    && make \
    && sudo make install
#~ PATH="/usr/local/bin/samtools/:$PATH"
########################################################################################
# cd to software dir
wget https://github.com/samtools/bcftools/releases/download/1.21/bcftools-1.21.tar.bz2 .
tar xjvf bcftools-1.21.tar.bz2
cd bcftools-1.21 \
    && ./configure \
    && make \
    && sudo make install
#~ PATH="/opt/bcftools/bin:$PATH"
########################################################################################
# cd to software dir
wget https://github.com/samtools/htslib/releases/download/1.21/htslib-1.21.tar.bz2 .
tar xjvf htslib-1.21.tar.bz2
cd htslib-1.21 \
    && ./configure \
    && make \
    && sudo make install
#~ PATH="/opt/htslib/bin:$PATH"
########################################################################################
# cd to software dir
wget https://github.com/freebayes/freebayes/releases/download/v1.3.6/freebayes-1.3.6-linux-amd64-static.gz .
gunzip freebayes-1.3.6-linux-amd64-static.gz \
    && chmod +x freebayes-1.3.6-linux-amd64-static \
    && sudo mkdir --parents /opt/freebayes/bin \
    && sudo mv freebayes-1.3.6-linux-amd64-static /opt/freebayes/bin/freebayes
#~ PATH="$PATH:/opt/freebayes/bin"
########################################################################################
wget https://github.com/lh3/minimap2/releases/download/v2.28/minimap2-2.28_x64-linux.tar.bz2 .
tar -jxvf minimap2-2.28_x64-linux.tar.bz2 \
    && sudo mv minimap2-2.28_x64-linux/ /opt/minimap2-2.28_x64-linux
#~ PATH="$PATH:/opt/minimap2-2.28_x64-linux"
########################################################################################
# Java:
#~ apt-get install -y default-jre
#~ JAVA_HOME=/usr/lib/jvm/default-java/
########################################################################################
# Canu: install from source...
sudo apt-get install -y libcurl4-gnutls-dev libpcap-dev libssl-dev
git clone https://github.com/marbl/canu.git \
    && cd canu/src \
    && make -j
#~ PATH="$PATH:/opt/canu/build/bin"
# ... or install latest canu release:
#~ wget https://github.com/marbl/canu/releases/download/v2.3/canu-2.3.Linux-amd64.tar.xz .
#~ tar -xJf canu-2.3.Linux-amd64.tar.xz \
    #~ && sudo mv canu-2.3/ /opt/canu-2.3
#~ PATH="/opt/canu-2.3/bin:$PATH"
# expected future patch release: https://github.com/marbl/canu/issues/2359#issuecomment-2631317846
########################################################################################
# For the above tools:
PATH="$PATH:/opt/freebayes/bin:/opt/minimap2-2.28_x64-linux:/opt/canu/build/bin"
########################################################################################
# For DNA Features Viewer:
pip install --no-cache-dir packaging
# Python packages:
pip install --no-cache-dir networkx==3.4.2 python-Levenshtein==0.26.1 pypugjs==5.11.0 jinja2==3.1.4 weasyprint==64.0 beautifulsoup4==4.12.3 Markdown==3.7
pip install --no-cache-dir matplotlib==3.10.0 numpy==2.2.0 pandas==2.2.3 weighted_levenshtein==0.2.2 portion==2.6.0
pip install --no-cache-dir biopython==1.84
pip install --no-cache-dir dna_features_viewer==3.1.3 geneblocks==1.2.3 cyvcf2==0.31.1 pdf_reports==0.3.7
pip install --no-cache-dir NanoPlot==1.44.1 nanofilt==2.8.0
pip install ediacara@git+https://github.com/Edinburgh-Genome-Foundry/Ediacara.git
pip install --no-cache-dir pypdf==5.3.0
pip install --no-cache-dir pydyf==0.11.0
