#!/usr/bin/env bash

SCRIPT_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )

# Go to script directory
cd ${SCRIPT_DIR}

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Download liftOver tool and the following chain files from UCSC:
# hg19ToMm10, mm9ToMm10, mm10Tohg19, and mm10Tohg38
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

if [ ! -f ${SCRIPT_DIR}/liftOver ]; then
    URL=https://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/liftOver
    curl -O ${URL}
    chmod u+x liftOver
fi

CHAIN="hg19ToMm10.over.chain"
if [ ! -f ${SCRIPT_DIR}/${CHAIN} ]; then
    URL=https://hgdownload.cse.ucsc.edu/goldenpath/hg19/liftOver/${CHAIN}.gz
    curl -O ${URL}
    gunzip ${CHAIN}.gz
fi

CHAIN="mm9ToMm10.over.chain"
if [ ! -f ${SCRIPT_DIR}/${CHAIN} ]; then
    URL=https://hgdownload.cse.ucsc.edu/goldenpath/mm9/liftOver/${CHAIN}.gz
    curl -O ${URL}
    gunzip ${CHAIN}.gz
fi


CHAIN="mm10ToHg19.over.chain"
if [ ! -f ${SCRIPT_DIR}/${CHAIN} ]; then
    URL=https://hgdownload.cse.ucsc.edu/goldenpath/mm10/liftOver/${CHAIN}.gz
    curl -O ${URL}
    gunzip ${CHAIN}.gz
fi


CHAIN="mm10ToHg38.over.chain"
if [ ! -f ${SCRIPT_DIR}/${CHAIN} ]; then
    URL=https://hgdownload.cse.ucsc.edu/goldenpath/mm10/liftOver/${CHAIN}.gz
    curl -O ${URL}
    gunzip ${CHAIN}.gz
fi

# Go back
cd - > /dev/null