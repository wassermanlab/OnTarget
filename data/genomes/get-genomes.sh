#!/usr/bin/env bash

SCRIPT_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )

# Go to script directory
cd ${SCRIPT_DIR}

for GENOME in "hg19" "hg38" "mm10"; do
    if [ ! -f ${GENOME}/${GENOME}.fa.sizes ]; then
        genomepy install -p UCSC -g ./ -t 8 -f ${GENOME}
    fi
done

# Go back
cd - > /dev/null