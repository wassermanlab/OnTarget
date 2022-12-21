#!/usr/bin/env bash

PY=/space/www/miniconda3/envs/ontarget/bin/genomepy

for GENOME in "hg19" "hg38" "mm10"; do
    if [ ! -f ${GENOME}/${GENOME}.fa.sizes ]; then
        ${PY} install -p UCSC -g ./ -t 8 -f ${GENOME}
    fi
done
