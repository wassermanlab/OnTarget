#!/usr/bin/env bash

# Get script directory
SCRIPT_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )

# Go to that directory
cd ${SCRIPT_DIR}

# Install several genomes using genomepy
for GENOME in "hg19" "hg38" "mm10"; do
    if [ ! -f ${GENOME}/${GENOME}.fa.sizes ]; then
        genomepy install -p UCSC -g ./ -t 8 -f ${GENOME}
    fi
done

# Go back
cd - > /dev/null