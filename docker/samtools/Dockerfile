FROM ubuntu:18.04

RUN apt-get update && apt-get install -y \
    wget \
    libbz2-dev \
    liblzma-dev \
    gcc \
    make \
    libncurses-dev \
    zlib1g-dev

WORKDIR /app

ENV HTSLIB_VER=1.10.2
ENV BCFTOOLS_VER=1.10.2
ENV SAMTOOLS_VER=1.10

# HTSLIB
RUN wget "https://github.com/samtools/htslib/releases/download/$HTSLIB_VER/htslib-$HTSLIB_VER.tar.bz2"
RUN tar -vxjf htslib-$HTSLIB_VER.tar.bz2
RUN cd htslib-$HTSLIB_VER &&\
    ./configure && \
    make && \
    make install
ENV PATH="$HOME/app/htslib-$HTSLIB_VER:$PATH"

# BCFTOOLS
RUN wget "https://github.com/samtools/bcftools/releases/download/$BCFTOOLS_VER/bcftools-$BCFTOOLS_VER.tar.bz2"
RUN tar -vxjf bcftools-$BCFTOOLS_VER.tar.bz2
RUN cd bcftools-$BCFTOOLS_VER && \
    ./configure && \
    make && \
    make install
ENV PATH="$HOME/app/bcftools-$BCFTOOLS_VER:$PATH"

# SAMTOOLS
RUN wget "https://github.com/samtools/samtools/releases/download/$SAMTOOLS_VER/samtools-$SAMTOOLS_VER.tar.bz2"
RUN tar -vxjf samtools-$SAMTOOLS_VER.tar.bz2
RUN cd samtools-$SAMTOOLS_VER && \
    ./configure && \
    make && \
    make install
ENV PATH="$HOME/app/samtools-$SAMTOOLS_VER:$PATH"

RUN rm samtools-$SAMTOOLS_VER.tar.bz2 && \
    rm bcftools-$BCFTOOLS_VER.tar.bz2 && \
    rm htslib-$HTSLIB_VER.tar.bz2

CMD /bin/bash
