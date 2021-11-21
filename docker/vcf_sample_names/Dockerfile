FROM python:3.8.6-slim

ENV BCFTOOLS_VER=1.10.2

# Install required packages
RUN apt-get update && apt-get install -y \
    wget gcc make bzip2 \
    libbz2-dev zlib1g-dev liblzma-dev libcurl4-gnutls-dev \
    libncurses-dev

# BCFTOOLS
RUN wget "https://github.com/samtools/bcftools/releases/download/$BCFTOOLS_VER/bcftools-$BCFTOOLS_VER.tar.bz2"
RUN tar -vxjf bcftools-$BCFTOOLS_VER.tar.bz2
RUN cd bcftools-$BCFTOOLS_VER && \
    ./configure && \
    make && \
    make install
ENV PATH="$HOME/app/bcftools-$BCFTOOLS_VER:$PATH"

# Copy Python script into container
ADD ./get_vcf_sample_names.py /app/get_vcf_sample_names.py
WORKDIR /app
CMD ["/bin/bash"]