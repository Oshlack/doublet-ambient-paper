#!/bin/bash

# This is a bash script that runs all doublet detection tools
# always remember to pipe out and error logs to a file while running:
# sh run_scds_scDblFinder_scrublet_doubletdetection.sh &> run_scds_scDblFinder_scrublet_doubletdetection.log

# inputs based on order
COUNTS=$1
OUTPUT=$2

mkdir ${OUTPUT}

echo "##########################################################################"
echo "Start doublet calling"
echo "Inputs detected:"
echo "  OUTPUT = ${OUTPUT}"
echo "  COUNTS = ${COUNTS}"
echo "##########################################################################"


#need to go here for renv reasons
cd /Users/meiermichelle/Library/CloudStorage/OneDrive-PeterMac/OSHLACK_LAB/PhD/ambient_doublet_singlecell

echo "##########################################################################"
echo "Start scds"
echo "##########################################################################"

# SCDS
Rscript ./code/doublet_tools/scds.R ${OUTPUT} ${COUNTS}

echo "##########################################################################"
echo "Start scDblFinder"
echo "##########################################################################"

# SCDBLFINDER
Rscript ./code/doublet_tools/scDblFinder.R ${OUTPUT} ${COUNTS}

# for mamba reasons
source /Users/meiermichelle/miniforge3/bin/activate
eval "$(mamba shell hook --shell bash)"

echo "##########################################################################"
echo "Start scrublet"
echo "##########################################################################"

# SCRUBLET
mamba activate scrublet
python ./code/doublet_tools/scrublet_tool.py ${OUTPUT} ${COUNTS}
mamba deactivate

echo "##########################################################################"
echo "Start DoubletDetection"
echo "##########################################################################"

# DOUBLETDETECTION
mamba activate DoubletDetection
python ./code/doublet_tools/DoubletDetection_tool.py ${OUTPUT} ${COUNTS}
mamba deactivate

echo "##########################################################################"
echo "Combine all doublet calling results"
echo "##########################################################################"

# COMBINE
Rscript ./code/doublet_tools/combine_results.R ${OUTPUT}

echo "##########################################################################"
echo "Doublet calling completed"
echo "##########################################################################"


