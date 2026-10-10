parent_dir="/Users/meiermichelle/Library/CloudStorage/OneDrive-PeterMac/OSHLACK_LAB/PhD/ambient_doublet_singlecell/output/data/Tonsils_10xFlex/Samples/"

cd /Users/meiermichelle/Library/CloudStorage/OneDrive-PeterMac/OSHLACK_LAB/PhD/ambient_doublet_singlecell

find "$parent_dir" -mindepth 1 -maxdepth 1 -type d | while read -r COUNTS; do

    sameName=$(basename "$COUNTS")

    echo "##########################################################################"
    echo "Start Sample ${sameName}"
    echo "##########################################################################"


    OUTPUT="${COUNTS}/doublet_calls" #make folder within

    mkdir ${OUTPUT}

    # Run
    sh ./code/doublet_tools/run_scds_scDblFinder_scrublet_doubletdetection.sh ${COUNTS} ${OUTPUT}

done
