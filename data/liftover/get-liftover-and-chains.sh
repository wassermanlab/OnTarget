#!/usr/bin/env bash

SCRIPT_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Download liftOver tool and the following chain files from UCSC:
# hg19ToMm10, mm9ToMm10, mm10Tohg19, and mm10Tohg38
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

if [ ! -f liftOver ]; then
    curl -o ${SCRIPT_DIR} \
        -O https://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/liftOver
    chmod u+x ${SCRIPT_DIR}/liftOver
fi

CHAIN="hg19ToMm10.over.chain"
if [ ! -f ${CHAIN} ]; then
    curl -o ${SCRIPT_DIR} \
        -O https://hgdownload.cse.ucsc.edu/goldenpath/hg19/liftOver/${CHAIN}.gz
    gunzip ${CHAIN}.gz
fi

CHAIN="mm9ToMm10.over.chain"
if [ ! -f ${CHAIN} ]; then
    curl -o ${SCRIPT_DIR} \
        -O https://hgdownload.cse.ucsc.edu/goldenpath/mm9/liftOver/${CHAIN}.gz
    gunzip ${CHAIN}.gz
fi

CHAIN="mm10ToHg19.over.chain"
if [ ! -f ${CHAIN} ]; then
    curl -o ${SCRIPT_DIR} \
        -O https://hgdownload.cse.ucsc.edu/goldenpath/mm10/liftOver/${CHAIN}.gz
    gunzip ${CHAIN}.gz
fi

CHAIN="mm10ToHg38.over.chain"
if [ ! -f ${CHAIN} ]; then
    curl -o ${SCRIPT_DIR} \
        -O https://hgdownload.cse.ucsc.edu/goldenpath/mm10/liftOver/${CHAIN}.gz
    gunzip ${CHAIN}.gz
fi
