#!/usr/bin/env bash

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Postnatal (day 0) mouse (mm10) midbrain processed DNase-seq peaks
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

DIR="dna_accessibility"

URL="https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM2195nnn"

GSM="GSM2195835"
ENC="ENCFF588FLL"
if [ ! -f ./${DIR}/${GSM}_${ENC}_peaks_mm10.bed.gz ]; then
    wget -P ${DIR} ${URL}/${GSM}/suppl/${GSM}_${ENC}_peaks_mm10.bed.gz
fi

GSM="GSM2195836"
ENC="ENCFF741TQE"
if [ ! -f ./${DIR}/${GSM}_${ENC}_peaks_mm10.bed.gz ]; then
    wget -P ${DIR} ${URL}/${GSM}/suppl/${GSM}_${ENC}_peaks_mm10.bed.gz
fi

if [ ! -f ./${DIR}.bed.gz ]; then
    zless ${DIR}/*.bed.gz | sort -k1,1 -k2,2n | \
        # awk '{print $1"\t"$2"\t"$3"\t"$1":"$2"-"$3;}' | \
        # mergeBed -c 4 -o collapse -delim "|" | gzip - > ./${DIR}.bed.gz
        awk '{print $1"\t"$2"\t"$3"\t"$1":"$2"-"$3;}' > ./${DIR}.bed.gz
fi

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Postnatal (day 0) mouse (mm10) midbrain processed histone ChIP-seq peaks
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

DIR="histone_modification"

URL="https://ftp.ncbi.nlm.nih.gov/geo/series/GSE82nnn"

GSM="GSE82687"
ENC="ENCFF770SEL"
HIST_MARK="H3K4me3"
if [ ! -f ./${DIR}/${HIST_MARK}/${GSM}_${ENC}_peaks_mm10.bed.gz ]; then
    wget -P ${DIR}/${HIST_MARK} \
        ${URL}/${GSM}/suppl/${GSM}_${ENC}_peaks_mm10.bed.gz
fi

if [ ! -f ./${DIR}-${HIST_MARK}.bed.gz ]; then
    zless ${DIR}/${HIST_MARK}/*.bed.gz | sort -k1,1 -k2,2n | \
        awk '{print $1"\t"$2"\t"$3"\t"$1":"$2"-"$3;}' | \
        gzip - > ./${DIR}-${HIST_MARK}.bed.gz
fi

GSM="GSE82870"
ENC="ENCFF991JUZ"
HIST_MARK="H3K27ac"
if [ ! -f ./${DIR}/${HIST_MARK}/${GSM}_${ENC}_peaks_mm10.bed.gz ]; then
    wget -P ${DIR}/${HIST_MARK} \
        ${URL}/${GSM}/suppl/${GSM}_${ENC}_peaks_mm10.bed.gz
fi

if [ ! -f ./${DIR}-${HIST_MARK}.bed.gz ]; then
    zless ${DIR}/${HIST_MARK}/*.bed.gz | sort -k1,1 -k2,2n | \
        awk '{print $1"\t"$2"\t"$3"\t"$1":"$2"-"$3;}' | \
        gzip - > ./${DIR}-${HIST_MARK}.bed.gz
fi

GSM="GSE82656"
ENC="ENCFF903IBK"
HIST_MARK="H3K4me1"
if [ ! -f ./${DIR}/${HIST_MARK}/${GSM}_${ENC}_peaks_mm10.bed.gz ]; then
    wget -P ${DIR}/${HIST_MARK} \
        ${URL}/${GSM}/suppl/${GSM}_${ENC}_peaks_mm10.bed.gz
fi

if [ ! -f ./${DIR}-${HIST_MARK}.bed.gz ]; then
    zless ${DIR}/${HIST_MARK}/*.bed.gz | sort -k1,1 -k2,2n | \
        awk '{print $1"\t"$2"\t"$3"\t"$1":"$2"-"$3;}' | \
        gzip - > ./${DIR}-${HIST_MARK}.bed.gz
fi

GSM="GSE82365"
ENC="ENCFF818XLK"
HIST_MARK="H3K9ac"
if [ ! -f ./${DIR}/${HIST_MARK}/${GSM}_${ENC}_peaks_mm10.bed.gz ]; then
    wget -P ${DIR}/${HIST_MARK} \
        ${URL}/${GSM}/suppl/${GSM}_${ENC}_peaks_mm10.bed.gz
fi

if [ ! -f ./${DIR}-${HIST_MARK}.bed.gz ]; then
    zless ${DIR}/${HIST_MARK}/*.bed.gz | sort -k1,1 -k2,2n | \
        awk '{print $1"\t"$2"\t"$3"\t"$1":"$2"-"$3;}' | \
        gzip - > ./${DIR}-${HIST_MARK}.bed.gz
fi

URL="https://ftp.ncbi.nlm.nih.gov/geo/series/GSE83nnn"

GSM="GSE83042"
ENC="ENCFF623CCW"
HIST_MARK="H3K36me3"
if [ ! -f ./${DIR}/${HIST_MARK}/${GSM}_${ENC}_peaks_mm10.bed.gz ]; then
    wget -P ${DIR}/${HIST_MARK} \
        ${URL}/${GSM}/suppl/${GSM}_${ENC}_peaks_mm10.bed.gz
fi

if [ ! -f ./${DIR}-${HIST_MARK}.bed.gz ]; then
    zless ${DIR}/${HIST_MARK}/*.bed.gz | sort -k1,1 -k2,2n | \
        awk '{print $1"\t"$2"\t"$3"\t"$1":"$2"-"$3;}' | \
        gzip - > ./${DIR}-${HIST_MARK}.bed.gz
fi
