#!/bin/bash

INPUT_DIR="/oshlack_lab/michelle.meier/projects/ambient_doublet_singlecell/output/data/Tonsils_10xFlex/Engineered/"

for INPUT_H5 in "${INPUT_DIR}"/Tonsils_10xFlex_engineered_12k_0.15p_*s_cellbenderInput_feature_bc_matrix.h5
do
    SAMPLE=$(basename "${INPUT_H5}" \
        "_cellbenderInput_feature_bc_matrix.h5")

    sbatch \
        --job-name="run_cellbender_${SAMPLE}" \
        --output="/oshlack_lab/michelle.meier/projects/ambient_doublet_singlecell/code/Tonsil_10xFlex/Engineered/cellbender_logs/run_CellBender_${SAMPLE}_%j.out" \
        --error="/oshlack_lab/michelle.meier/projects/ambient_doublet_singlecell/code/Tonsil_10xFlex/Engineered/cellbender_logs/run_CellBender_${SAMPLE}_%j.err" \
        run_CellBender_EngineeredDatasets.sbatch "${INPUT_H5}" "${SAMPLE}"

done
