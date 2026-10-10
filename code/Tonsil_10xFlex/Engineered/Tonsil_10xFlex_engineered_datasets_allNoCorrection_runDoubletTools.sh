#!/bin/bash

INPUT="/Users/meiermichelle/Library/CloudStorage/OneDrive-PeterMac/OSHLACK_LAB/PhD/ambient_doublet_singlecell/output/data/Tonsils_10xFlex/Engineered/"

cd /Users/meiermichelle/Library/CloudStorage/OneDrive-PeterMac/OSHLACK_LAB/PhD/ambient_doublet_singlecell

find "$INPUT" -mindepth 1 -maxdepth 1 -type d -name "Tonsils_10xFlex_engineered_12k_0.15p*s" | while read -r COUNTS; do

    sameName=$(basename "$COUNTS")

    echo "##########################################################################"
    echo "Start Sample ${sameName}"
    echo "##########################################################################"


    OUTPUT="${COUNTS}/doublet_calls" #make folder within

    mkdir ${OUTPUT}

    # Run
    sh ./code/doublet_tools/run_cxds_bcds_scds_scDblFinder_scrublet_doubletdetection.sh ${COUNTS} ${OUTPUT}

done
