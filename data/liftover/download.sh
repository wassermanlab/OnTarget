#!/usr/bin/env bash

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Download liftOver and chains (hg19ToMm10, mm9ToMm10, mm10Tohg19, and
# mm10Tohg38) from UCSC
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

if [ ! -f liftOver ]; then
    wget http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/liftOver
    chmod u+x liftOver
fi

CHAIN="hg19ToMm10.over.chain.gz"
if [ ! -f ${CHAIN} ]; then
    wget https://hgdownload.cse.ucsc.edu/goldenpath/hg19/liftOver/${CHAIN}
fi

CHAIN="mm9ToMm10.over.chain.gz"
if [ ! -f ${CHAIN} ]; then
    wget https://hgdownload.cse.ucsc.edu/goldenpath/mm9/liftOver/${CHAIN}
fi

CHAIN="mm10ToHg19.over.chain.gz"
if [ ! -f ${CHAIN} ]; then
    wget https://hgdownload.cse.ucsc.edu/goldenpath/mm10/liftOver/${CHAIN}
fi

CHAIN="mm10ToHg38.over.chain.gz"
if [ ! -f ${CHAIN} ]; then
    wget https://hgdownload.cse.ucsc.edu/goldenpath/mm10/liftOver/${CHAIN}
fi