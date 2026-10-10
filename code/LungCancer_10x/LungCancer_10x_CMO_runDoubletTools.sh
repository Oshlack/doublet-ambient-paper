#!/bin/bash

INPUT="/Users/meiermichelle/Library/CloudStorage/OneDrive-PeterMac/OSHLACK_LAB/PhD/ambient_doublet_singlecell/output/data/LungCancer_10x/"

cd /Users/meiermichelle/Library/CloudStorage/OneDrive-PeterMac/OSHLACK_LAB/PhD/ambient_doublet_singlecell

find "$INPUT" -mindepth 1 -maxdepth 1 -type d -name "LungCancer_10x_cellranger_allmultip*" | while read -r COUNTS; do

    sameName=$(basename "$COUNTS")

    echo "##########################################################################"
    echo "Start Dataset ${sameName}"
    echo "##########################################################################"


    OUTPUT="${COUNTS}/doublet_calls" #make folder within

    mkdir ${OUTPUT}

    # Run
    sh ./code/doublet_tools/run_cxds_bcds_scds_scDblFinder_scrublet_doubletdetection_filtered.sh ${COUNTS} ${OUTPUT}

done
