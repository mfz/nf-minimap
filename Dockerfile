FROM --platform=linux/amd64 ubuntu:22.04

ENV DEBIAN_FRONTEND=noninteractive
ENV SAMTOOLS_VERSION=1.23
ENV HTSLIB_VERSION=1.23
ENV BAMTOOLS_VERSION=2.5.2
ENV MINIMAP_VERSION=2.28

RUN sed -i 's|http://archive.ubuntu.com/ubuntu/|http://ubuntu.hysing.is/ubuntu/|g' /etc/apt/sources.list \
    && sed -i 's|http://security.ubuntu.com/ubuntu/|http://ubuntu.hysing.is/ubuntu/|g' /etc/apt/sources.list

# Install build dependencies (pinned via Ubuntu release)
RUN apt-get update && apt-get install -y \
    build-essential=12.9ubuntu3 \
    wget \
    curl \
    ca-certificates \
    git \
    python3 \
    python3-pip \
    python-is-python3 \
    autoconf \
    automake \
    pkg-config \
    cmake \
    zlib1g-dev \
    libbz2-dev \
    liblzma-dev \
    libcurl4-openssl-dev \
    libssl-dev \
    libncurses5-dev \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /opt

COPY view_region.py /usr/local/bin/view_region.py
COPY process.sh /usr/local/bin/process.sh
RUN chmod 755 /usr/local/bin/view_region.py \
    && chmod 755 /usr/local/bin/process.sh \
    && pip3 install --no-cache-dir pysam

# Build htslib
RUN wget https://github.com/samtools/htslib/releases/download/${HTSLIB_VERSION}/htslib-${HTSLIB_VERSION}.tar.bz2 \
    && tar -xjf htslib-${HTSLIB_VERSION}.tar.bz2 \
    && cd htslib-${HTSLIB_VERSION} \
    && ./configure --enable-libcurl \
    && make -j4 \
    && make install

# Build samtools
RUN wget https://github.com/samtools/samtools/releases/download/${SAMTOOLS_VERSION}/samtools-${SAMTOOLS_VERSION}.tar.bz2 \
    && tar -xjf samtools-${SAMTOOLS_VERSION}.tar.bz2 \
    && cd samtools-${SAMTOOLS_VERSION} \
    && ./configure \
    && make -j4 \
    && make install


# Install minimap2    
RUN wget https://github.com/lh3/minimap2/releases/download/v${MINIMAP_VERSION}/minimap2-${MINIMAP_VERSION}_x64-linux.tar.bz2 \
    && tar -xjf minimap2-${MINIMAP_VERSION}_x64-linux.tar.bz2 \
    && cp minimap2-${MINIMAP_VERSION}_x64-linux/minimap2 /usr/local/bin/



# Build bamtools
RUN wget https://github.com/pezmaster31/bamtools/archive/refs/tags/v${BAMTOOLS_VERSION}.tar.gz \
    && tar -xzf v${BAMTOOLS_VERSION}.tar.gz \
    && cd bamtools-${BAMTOOLS_VERSION} \
    && mkdir build \
    && cd build \
    && cmake -DCMAKE_INSTALL_PREFIX=/usr/local .. \
    && make -j4 \
    && make install

# Build Telseq from pinned commit
RUN git clone https://github.com/zd1/telseq.git \
    && cd telseq/src \
    && ./autogen.sh \
    && ./configure --with-bamtools=/usr/local \
    && make \
    && cp Telseq/telseq /usr/local/bin/

WORKDIR /data

CMD ["/bin/bash"]
