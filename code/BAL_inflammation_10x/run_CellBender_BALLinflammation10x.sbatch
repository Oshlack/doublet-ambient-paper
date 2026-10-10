#!/bin/bash
#SBATCH --partition=rhel_gpu
#SBATCH --gres=gpu:T4:1
#SBATCH --cpus-per-task=1
#SBATCH --ntasks=1
#SBATCH --mem=48G
#SBATCH --time=00-06:00:00
#SBATCH --job-name="run_cellbender_BAL_10x"
#SBATCH --output="/oshlack_lab/michelle.meier/projects/ambient_doublet_singlecell/code/BAL_inflammation_10x/cellbender_logs/run_cellbender_BAL_10x_cellranger_%j.out" \
#SBATCH --error="/oshlack_lab/michelle.meier/projects/ambient_doublet_singlecell/code/BAL_inflammation_10x/cellbender_logs/run_cellbender_BAL_10x_cellranger_%j.err" \

set -euo pipefail

# INPUT
INPUT_H5="/oshlack_lab/michelle.meier/projects/ambient_doublet_singlecell/data/BAL_inflammation_10x/C133_batch3_1/raw_feature_bc_matrix.h5"
SAMPLE="BAL_inflammation_10x"
SAMPLE_OUT="${SAMPLE}_CellBender"

# Paths
OUTPUT_DIR="/oshlack_lab/michelle.meier/projects/ambient_doublet_singlecell/output/data/BAL_inflammation_10x/${SAMPLE_OUT}"

# make output
mkdir -p "${OUTPUT_DIR}"

# Change to output directory
cd "${OUTPUT_DIR}"

# Load required modules
module load cuda/12.4.0-gcc-13.2.0   # for RHEL9
module load miniconda3/24.3.0-gcc-13.2.0   # for RHEL9

# Activate CellBender conda environment
source /config/spack/v0.22.1/opt/spack/linux-rhel9-haswell/gcc-13.2.0/miniconda3-24.3.0-re4lxkwqjaazi6ckedtfhnzrjiefmhei/etc/profile.d/conda.sh
conda activate cellbender_v3.0.2

# Run CellBender remove‑background
cellbender remove-background \
    --input  "${INPUT_H5}" \
    --output "${OUTPUT_DIR}/${SAMPLE_OUT}.h5" \
    --expected-cells 33000 \
    --total-droplets-included 65000 \
    --exclude-feature-types "Multiplexing Capture" \
    --cuda

# Make compatible with seurat
ptrepack --complevel 5 ${SAMPLE_OUT}.h5:/matrix ${SAMPLE_OUT}_seurat.h5:/matrix

echo "✅ CellBender completed for ${SAMPLE}. Outputs in ${OUTPUT_DIR}"
