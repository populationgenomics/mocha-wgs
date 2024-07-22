FROM debian:bookworm-slim

ENV PATH /software/htslib-1.20/build/bin:/software/bcftools-1.20/build/bin:$PATH

RUN apt-get update && apt-get install -y libbz2-dev liblzma-dev git wget bash bzip2 zip curl gcc zlib1g-dev gawk build-essential \
  && rm -r /var/lib/apt/lists/* \
  && rm -r /var/cache/apt/* \
  && mkdir -p /software/tmp \
  && mkdir -p /software/htslib-1.20 \
  && mkdir -p /software/bcftools-1.20 \
  && cd /software/tmp \
  && wget https://github.com/samtools/htslib/releases/download/1.20/htslib-1.20.tar.bz2 \
  && wget https://github.com/samtools/bcftools/releases/download/1.20/bcftools-1.20.tar.bz2 \
  && tar -xjvf htslib-1.20.tar.bz2 \
  && tar -xjvf bcftools-1.20.tar.bz2 \
  && cd /software/tmp/htslib-1.20 \
  && ./configure --prefix=/software/htslib-1.20/build \
  && make \
  && make install \
  && cd - \
  && cd /software/tmp/bcftools-1.20 \
  && for f in mocha.h beta_binom.h genome_rules.h mocha.c trio-phase.c mochatools.c extendFMT.c; do wget -P plugins https://raw.githubusercontent.com/freeseek/mocha/master/${f}; done \
  && ./configure --prefix=/software/bcftools-1.20/build \
  && make \
  && make install \
  && cd /software \
  && rm -rf /software/tmp
