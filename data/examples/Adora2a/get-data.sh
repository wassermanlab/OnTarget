#!/usr/bin/env bash

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Postmortem human (hg19) striatal neurons (both from the accumbens nucleus
# and putamen) processed ATAC-seq data (BOCA; PMID: XXXXX); liftOver mm10
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

DIR="dna_accessibility"

URL="https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM2546nnn"

GSM="GSM2546463"
BOCA="372_N_NAC"

if [ ! -f ./${DIR}/${GSM}_${BOCA}_norm_sorted_rmdup_nochrM.bw ]; then
    wget -P ${DIR} ${URL}/${GSM}/suppl/${GSM}_${BOCA}_norm_sorted_rmdup_nochrM.bw
fi

if [ ! -f ./${DIR}/${GSM}_${BOCA}.bg ]; then
    bigWigToBedGraph ./${DIR}/${GSM}_${BOCA}_norm_sorted_rmdup_nochrM.bw \
        ./${DIR}/${GSM}_${BOCA}.bg
fi

if [ ! -f ./${DIR}/${GSM}_${BOCA}.bed ]; then
    macs2 bdgpeakcall -i ./${DIR}/${GSM}_${BOCA}.bg \
        -o ./${DIR}/${GSM}_${BOCA}.bed 
fi

GSM="GSM2546440"
BOCA="351_N_PUT"
if [ ! -f ./${DIR}/${GSM}_${BOCA}_norm_sorted_rmdup_nochrM.bw ]; then
    wget -P ${DIR} ${URL}/${GSM}/suppl/${GSM}_${BOCA}_norm_sorted_rmdup_nochrM.bw
fi

if [ ! -f ./${DIR}/${GSM}_${BOCA}.bg ]; then
    bigWigToBedGraph ./${DIR}/${GSM}_${BOCA}_norm_sorted_rmdup_nochrM.bw \
        ./${DIR}/${GSM}_${BOCA}.bg
fi

if [ ! -f ./${DIR}/${GSM}_${BOCA}.bed ]; then
    macs2 bdgpeakcall -i ./${DIR}/${GSM}_${BOCA}.bg \
        -o ./${DIR}/${GSM}_${BOCA}.bed 
fi

GSM="GSM2546465"
BOCA="372_N_PUT"
if [ ! -f ./${DIR}/${GSM}_${BOCA}_norm_sorted_rmdup_nochrM.bw ]; then
    wget -P ${DIR} ${URL}/${GSM}/suppl/${GSM}_${BOCA}_norm_sorted_rmdup_nochrM.bw
fi

if [ ! -f ./${DIR}/${GSM}_${BOCA}.bg ]; then
    bigWigToBedGraph ./${DIR}/${GSM}_${BOCA}_norm_sorted_rmdup_nochrM.bw \
        ./${DIR}/${GSM}_${BOCA}.bg
fi

if [ ! -f ./${DIR}/${GSM}_${BOCA}.bed ]; then
    macs2 bdgpeakcall -i ./${DIR}/${GSM}_${BOCA}.bg \
        -o ./${DIR}/${GSM}_${BOCA}.bed 
fi

GSM="GSM2546489"
BOCA="437_N_PUT"
if [ ! -f ./${DIR}/${GSM}_${BOCA}_norm_sorted_rmdup_nochrM.bw ]; then
    wget -P ${DIR} ${URL}/${GSM}/suppl/${GSM}_${BOCA}_norm_sorted_rmdup_nochrM.bw
fi

if [ ! -f ./${DIR}/${GSM}_${BOCA}.bg ]; then
    bigWigToBedGraph ./${DIR}/${GSM}_${BOCA}_norm_sorted_rmdup_nochrM.bw \
        ./${DIR}/${GSM}_${BOCA}.bg
fi

if [ ! -f ./${DIR}/${GSM}_${BOCA}.bed ]; then
    macs2 bdgpeakcall -i ./${DIR}/${GSM}_${BOCA}.bg \
        -o ./${DIR}/${GSM}_${BOCA}.bed 
fi

GSM="GSM2546535"
BOCA="446_N_PUT"
if [ ! -f ./${DIR}/${GSM}_${BOCA}_norm_sorted_rmdup_nochrM.bw ]; then
    wget -P ${DIR} ${URL}/${GSM}/suppl/${GSM}_${BOCA}_norm_sorted_rmdup_nochrM.bw
fi

if [ ! -f ./${DIR}/${GSM}_${BOCA}.bg ]; then
    bigWigToBedGraph ./${DIR}/${GSM}_${BOCA}_norm_sorted_rmdup_nochrM.bw \
        ./${DIR}/${GSM}_${BOCA}.bg
fi

if [ ! -f ./${DIR}/${GSM}_${BOCA}.bed ]; then
    macs2 bdgpeakcall -i ./${DIR}/${GSM}_${BOCA}.bg \
        -o ./${DIR}/${GSM}_${BOCA}.bed 
fi

if [ ! -f ./${DIR}.bed.gz ]; then
    less ${DIR}/GSM2546*.bed | sort -k1,1 -k2,2n | \
        awk '{print $1"\t"$2"\t"$3"\t"$1":"$2"-"$3;}' > ${DIR}/hg19.bed
    ../../chains/liftOver -minMatch=0.1 ${DIR}/hg19.bed \
        ../../chains/hg19ToMm10.over.chain.gz \
        ./${DIR}.bed ${DIR}/unMapped
    gzip ./${DIR}.bed
fi

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#  Adult (8 weeks old) mouse (mm9) accumbens nucleus processed
# (saline treatment only) histone (H3K4me1, H3K4me3, H3K36me3) and Pol2ra
# ChIP-seq data in bigWig format
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

DIR="histone_modification"

URL="https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM1050nnn"

HIST_MARK="H3K4me1"

GSM="GSM1050349"
if [ ! -f ./${DIR}/${HIST_MARK}/${GSM}_mm9.nac.h3k4me1.sal.rep1.bw ]; then
    wget -P ${DIR}/${HIST_MARK} \
        ${URL}/${GSM}/suppl/${GSM}_mm9.nac.h3k4me1.sal.rep1.bw
fi

if [ ! -f ./${DIR}/${HIST_MARK}/rep1.bg ]; then
    bigWigToBedGraph \
        ./${DIR}/${HIST_MARK}/${GSM}_mm9.nac.h3k4me1.sal.rep1.bw \
        ./${DIR}/${HIST_MARK}/rep1.bg
fi

if [ ! -f ./${DIR}/${HIST_MARK}/rep1.bed ]; then
    macs2 bdgpeakcall -i ./${DIR}/${HIST_MARK}/rep1.bg \
        -o ./${DIR}/${HIST_MARK}/rep1.bed
fi

GSM="GSM1050350"
if [ ! -f ./${DIR}/${HIST_MARK}/${GSM}_mm9.nac.h3k4me1.sal.rep2.bw ]; then
    wget -P ${DIR}/${HIST_MARK} \
        ${URL}/${GSM}/suppl/${GSM}_mm9.nac.h3k4me1.sal.rep2.bw
fi

if [ ! -f ./${DIR}/${HIST_MARK}/rep2.bg ]; then
    bigWigToBedGraph \
        ./${DIR}/${HIST_MARK}/${GSM}_mm9.nac.h3k4me1.sal.rep2.bw \
        ./${DIR}/${HIST_MARK}/rep2.bg
fi

if [ ! -f ./${DIR}/${HIST_MARK}/rep2.bed ]; then
    macs2 bdgpeakcall -i ./${DIR}/${HIST_MARK}/rep2.bg \
        -o ./${DIR}/${HIST_MARK}/rep2.bed
fi

GSM="GSM1050351"
if [ ! -f ./${DIR}/${HIST_MARK}/${GSM}_mm9.nac.h3k4me1.sal.rep3.bw ]; then
    wget -P ${DIR}/${HIST_MARK} \
        ${URL}/${GSM}/suppl/${GSM}_mm9.nac.h3k4me1.sal.rep3.bw
fi

if [ ! -f ./${DIR}/${HIST_MARK}/rep3.bg ]; then
    bigWigToBedGraph \
         ./${DIR}/${HIST_MARK}/${GSM}_mm9.nac.h3k4me1.sal.rep3.bw \
        ./${DIR}/${HIST_MARK}/rep3.bg
fi

if [ ! -f ./${DIR}/${HIST_MARK}/rep3.bed ]; then
    macs2 bdgpeakcall -i ./${DIR}/${HIST_MARK}/rep3.bg \
        -o ./${DIR}/${HIST_MARK}/rep3.bed
fi

if [ ! -f ./${DIR}-${HIST_MARK}.bed.gz ]; then
    less ./${DIR}/${HIST_MARK}/rep*.bed | sort -k1,1 -k2,2n | \
        awk '{print $1"\t"$2"\t"$3"\t"$1":"$2"-"$3;}' \
        > ./${DIR}/${HIST_MARK}/mm9.bed
    ../../chains/liftOver ./${DIR}/${HIST_MARK}/mm9.bed \
        ../../chains/mm9ToMm10.over.chain.gz \
        ./${DIR}-${HIST_MARK}.bed ./${DIR}/${HIST_MARK}/unMapped
    gzip ./${DIR}-${HIST_MARK}.bed
fi

HIST_MARK="H3K4me3"

GSM="GSM1050355"
if [ ! -f ./${DIR}/${HIST_MARK}/${GSM}_mm9.nac.h3k4me3.sal.rep1.bw ]; then
    wget -P ${DIR}/${HIST_MARK} \
        ${URL}/${GSM}/suppl/${GSM}_mm9.nac.h3k4me3.sal.rep1.bw
fi

if [ ! -f ./${DIR}/${HIST_MARK}/rep1.bg ]; then
    bigWigToBedGraph \
        ./${DIR}/${HIST_MARK}/${GSM}_mm9.nac.h3k4me3.sal.rep1.bw \
        ./${DIR}/${HIST_MARK}/rep1.bg
fi

if [ ! -f ./${DIR}/${HIST_MARK}/rep1.bed ]; then
    macs2 bdgpeakcall -i ./${DIR}/${HIST_MARK}/rep1.bg \
        -o ./${DIR}/${HIST_MARK}/rep1.bed
fi

GSM="GSM1050356"
if [ ! -f ./${DIR}/${HIST_MARK}/${GSM}_mm9.nac.h3k4me3.sal.rep2.bw ]; then
    wget -P ${DIR}/${HIST_MARK} \
        ${URL}/${GSM}/suppl/${GSM}_mm9.nac.h3k4me3.sal.rep2.bw
fi

if [ ! -f ./${DIR}/${HIST_MARK}/rep2.bg ]; then
    bigWigToBedGraph \
        ./${DIR}/${HIST_MARK}/${GSM}_mm9.nac.h3k4me3.sal.rep2.bw \
        ./${DIR}/${HIST_MARK}/rep2.bg
fi

if [ ! -f ./${DIR}/${HIST_MARK}/rep2.bed ]; then
    macs2 bdgpeakcall -i ./${DIR}/${HIST_MARK}/rep2.bg \
        -o ./${DIR}/${HIST_MARK}/rep2.bed
fi

GSM="GSM1050357"
if [ ! -f ./${DIR}/${HIST_MARK}/${GSM}_mm9.nac.h3k4me3.sal.rep3.bw ]; then
    wget -P ${DIR}/${HIST_MARK} \
        ${URL}/${GSM}/suppl/${GSM}_mm9.nac.h3k4me3.sal.rep3.bw
fi

if [ ! -f ./${DIR}/${HIST_MARK}/rep3.bg ]; then
    bigWigToBedGraph \
         ./${DIR}/${HIST_MARK}/${GSM}_mm9.nac.h3k4me3.sal.rep3.bw \
        ./${DIR}/${HIST_MARK}/rep3.bg
fi

if [ ! -f ./${DIR}/${HIST_MARK}/rep3.bed ]; then
    macs2 bdgpeakcall -i ./${DIR}/${HIST_MARK}/rep3.bg \
        -o ./${DIR}/${HIST_MARK}/rep3.bed
fi

if [ ! -f ./${DIR}-${HIST_MARK}.bed.gz ]; then
    less ./${DIR}/${HIST_MARK}/rep*.bed | sort -k1,1 -k2,2n | \
        awk '{print $1"\t"$2"\t"$3"\t"$1":"$2"-"$3;}' \
        > ./${DIR}/${HIST_MARK}/mm9.bed
    ../../chains/liftOver ./${DIR}/${HIST_MARK}/mm9.bed \
        ../../chains/mm9ToMm10.over.chain.gz \
        ./${DIR}-${HIST_MARK}.bed ./${DIR}/${HIST_MARK}/unMapped
    gzip ./${DIR}-${HIST_MARK}.bed
fi

HIST_MARK="H3K36me3"

GSM="GSM1050343"
if [ ! -f ./${DIR}/${HIST_MARK}/${GSM}_mm9.nac.h3k36me3.sal.rep1.bw ]; then
    wget -P ${DIR}/${HIST_MARK} \
        ${URL}/${GSM}/suppl/${GSM}_mm9.nac.h3k36me3.sal.rep1.bw
fi

if [ ! -f ./${DIR}/${HIST_MARK}/rep1.bg ]; then
    bigWigToBedGraph \
        ./${DIR}/${HIST_MARK}/${GSM}_mm9.nac.h3k36me3.sal.rep1.bw \
        ./${DIR}/${HIST_MARK}/rep1.bg
fi

if [ ! -f ./${DIR}/${HIST_MARK}/rep1.bed ]; then
    macs2 bdgpeakcall -i ./${DIR}/${HIST_MARK}/rep1.bg \
        -o ./${DIR}/${HIST_MARK}/rep1.bed
fi

GSM="GSM1050344"
if [ ! -f ./${DIR}/${HIST_MARK}/${GSM}_mm9.nac.h3k36me3.sal.rep2.bw ]; then
    wget -P ${DIR}/${HIST_MARK} \
        ${URL}/${GSM}/suppl/${GSM}_mm9.nac.h3k36me3.sal.rep2.bw
fi

if [ ! -f ./${DIR}/${HIST_MARK}/rep2.bg ]; then
    bigWigToBedGraph \
        ./${DIR}/${HIST_MARK}/${GSM}_mm9.nac.h3k36me3.sal.rep2.bw \
        ./${DIR}/${HIST_MARK}/rep2.bg
fi

if [ ! -f ./${DIR}/${HIST_MARK}/rep2.bed ]; then
    macs2 bdgpeakcall -i ./${DIR}/${HIST_MARK}/rep2.bg \
        -o ./${DIR}/${HIST_MARK}/rep2.bed
fi

GSM="GSM1050345"
if [ ! -f ./${DIR}/${HIST_MARK}/${GSM}_mm9.nac.h3k36me3.sal.rep3.bw ]; then
    wget -P ${DIR}/${HIST_MARK} \
        ${URL}/${GSM}/suppl/${GSM}_mm9.nac.h3k36me3.sal.rep3.bw
fi

if [ ! -f ./${DIR}/${HIST_MARK}/rep3.bg ]; then
    bigWigToBedGraph \
         ./${DIR}/${HIST_MARK}/${GSM}_mm9.nac.h3k36me3.sal.rep3.bw \
        ./${DIR}/${HIST_MARK}/rep3.bg
fi

if [ ! -f ./${DIR}/${HIST_MARK}/rep3.bed ]; then
    macs2 bdgpeakcall -i ./${DIR}/${HIST_MARK}/rep3.bg \
        -o ./${DIR}/${HIST_MARK}/rep3.bed
fi

if [ ! -f ./${DIR}-${HIST_MARK}.bed.gz ]; then
    less ./${DIR}/${HIST_MARK}/rep*.bed | sort -k1,1 -k2,2n | \
        awk '{print $1"\t"$2"\t"$3"\t"$1":"$2"-"$3;}' \
        > ./${DIR}/${HIST_MARK}/mm9.bed
    ../../chains/liftOver ./${DIR}/${HIST_MARK}/mm9.bed \
        ../../chains/mm9ToMm10.over.chain.gz \
        ./${DIR}-${HIST_MARK}.bed ./${DIR}/${HIST_MARK}/unMapped
    gzip ./${DIR}-${HIST_MARK}.bed
fi

DIR="tf_binding"

TF="POL2"

GSM="GSM1050368"
if [ ! -f ./${DIR}/${TF}/${GSM}_mm9.nac.pol2.sal.rep1.bw ]; then
    wget -P ${DIR}/${TF} \
        ${URL}/${GSM}/suppl/${GSM}_mm9.nac.pol2.sal.rep1.bw
fi

if [ ! -f ./${DIR}/${TF}/rep1.bg ]; then
    bigWigToBedGraph \
        ./${DIR}/${TF}/${GSM}_mm9.nac.pol2.sal.rep1.bw \
        ./${DIR}/${TF}/rep1.bg
fi

if [ ! -f ./${DIR}/${TF}/rep1.bed ]; then
    macs2 bdgpeakcall -i ./${DIR}/${TF}/rep1.bg \
        -o ./${DIR}/${TF}/rep1.bed
fi

GSM="GSM1050369"
if [ ! -f ./${DIR}/${TF}/${GSM}_mm9.nac.pol2.sal.rep2.bw ]; then
    wget -P ${DIR}/${TF} \
        ${URL}/${GSM}/suppl/${GSM}_mm9.nac.pol2.sal.rep2.bw
fi

if [ ! -f ./${DIR}/${TF}/rep2.bg ]; then
    bigWigToBedGraph \
        ./${DIR}/${TF}/${GSM}_mm9.nac.pol2.sal.rep2.bw \
        ./${DIR}/${TF}/rep2.bg
fi

if [ ! -f ./${DIR}/${TF}/rep2.bed ]; then
    macs2 bdgpeakcall -i ./${DIR}/${TF}/rep2.bg \
        -o ./${DIR}/${TF}/rep2.bed
fi

GSM="GSM1050370"
if [ ! -f ./${DIR}/${TF}/${GSM}_mm9.nac.pol2.sal.rep3.bw ]; then
    wget -P ${DIR}/${TF} \
        ${URL}/${GSM}/suppl/${GSM}_mm9.nac.pol2.sal.rep3.bw
fi

if [ ! -f ./${DIR}/${TF}/rep3.bg ]; then
    bigWigToBedGraph \
         ./${DIR}/${TF}/${GSM}_mm9.nac.pol2.sal.rep3.bw \
        ./${DIR}/${TF}/rep3.bg
fi

if [ ! -f ./${DIR}/${TF}/rep3.bed ]; then
    macs2 bdgpeakcall -i ./${DIR}/${TF}/rep3.bg \
        -o ./${DIR}/${TF}/rep3.bed
fi

if [ ! -f ./${DIR}-${TF}.bed.gz ]; then
    less ./${DIR}/${TF}/rep*.bed | sort -k1,1 -k2,2n | \
        awk '{print $1"\t"$2"\t"$3"\t"$1":"$2"-"$3;}' > ./${DIR}/${TF}/mm9.bed
    ../../chains/liftOver ./${DIR}/${HIST_MARK}/mm9.bed \
        ../../chains/mm9ToMm10.over.chain.gz \
        ./${DIR}-${HIST_MARK}.bed ./${DIR}/${HIST_MARK}/unMapped
    gzip ./${DIR}-${HIST_MARK}.bed
fi

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#  Mouse (mm9) medium spiny neurons processed H3K27ac ChIP-seq peaks
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

DIR="histone_modification"

URL="https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM2230nnn/GSM2230267/suppl"
BED_FILE="GSM2230267_H3K27ac_neuron_over_input.3kb_up_down_tss.bed.gz"
HIST_MARK="H3K27ac"
if [ ! -f ./${DIR}/${HIST_MARK}/${BED_FILE} ]; then
    wget -P ${DIR}/${HIST_MARK} ${URL}/${BED_FILE}
fi

if [ ! -f ./${DIR}-${HIST_MARK}.bed.gz ]; then
    zless ./${DIR}/${HIST_MARK}/${BED_FILE} | cut -f 1-3 | \
        sort -k1,1 -k2,2n | uniq > ./${DIR}/${HIST_MARK}/mm9.bed
    ../../chains/liftOver ./${DIR}/${HIST_MARK}/mm9.bed \
        ../../chains/mm9ToMm10.over.chain.gz \
        ./${DIR}-${HIST_MARK}.bed ./${DIR}/${HIST_MARK}/unMapped
    gzip ./${DIR}-${HIST_MARK}.bed
fi
