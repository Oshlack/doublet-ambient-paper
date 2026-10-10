#!/usr/env/bin python

################################################################################
# RUN SCRIPT FOR SCRUBLET DOUBLET CALLING
################################################################################
# reference: https://anaconda.org/channels/bioconda/packages/scrublet/overview
# this script is adapted from demuxafy: https://demultiplexing-doublet-detecting-docs.readthedocs.io/en/latest/Scrublet.html

# load library
import sys
import os
import scrublet as scr
import scipy.io
import matplotlib
matplotlib.use('AGG')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import umap
import numba
import numba.typed
import argparse

# Get path of mods directory from current script directory
mods_path = "/Users/meiermichelle/Library/CloudStorage/OneDrive-PeterMac/OSHLACK_LAB/PhD/ambient_doublet_singlecell/code/doublet_tools/" ## this is the path to the custom script for loading 10x data in python
sys.path.append(mods_path)
import read10x


## Set up parameters from trailing inputs
parser = argparse.ArgumentParser()

parser.add_argument("outdir", help="Output directory")
parser.add_argument("counts_matrix_dir", help="Path to counts matrix directory")


args = parser.parse_args()

counts_matrix_dir = args.counts_matrix_dir
outdir = args.outdir


if not os.path.isdir(outdir):
  os.mkdir(outdir)
  
plt.rc('font', size=14)
plt.rcParams['pdf.fonttype'] = 42

## Basic run with scrublet
counts_matrix = read10x.import_cellranger_mtx(counts_matrix_dir) 

try:
  barcodes_df = read10x.read_barcodes(counts_matrix_dir + "/barcodes.tsv.gz")
except:
  try:
    barcodes_df = read10x.read_barcodes(counts_matrix_dir + "/barcodes.tsv")
  except:
    print("No barcode file in provided counts matrix directory. Please double check the directory or provide the full path to the barcode file to use.")

#read filtered barcodes
barcodes_filtered_df = pd.read_csv(counts_matrix_dir + "/filtered_qc_barcodes.tsv", 
sep = "\t", header=1, names = ["Barcode"])

#subset + reorder to match filtered order
barcode_to_idx = pd.Series(
    np.arange(len(barcodes_df)),
    index=barcodes_df['Barcode']
)

row_idx = barcode_to_idx.loc[
    barcodes_filtered_df['Barcode']
].to_numpy()

counts_matrix = counts_matrix[row_idx, :]


#dbl_rate = counts_matrix.shape[0]/1000 * 0.008 ## This is the calculation for 10x doublet rate but will be different for other platforms
# use default doublet rate instead, which is 0.1 

print('Counts matrix shape: {} rows, {} columns'.format(counts_matrix.shape[0], counts_matrix.shape[1]))
scrub = scr.Scrublet(counts_matrix)
doublet_scores, predicted_doublets = scrub.scrub_doublets(min_counts=3,
                                                          min_cells=3,
                                                          min_gene_variability_pctl=85, #default, but may need to be changed if distributions look bad
                                                          n_prin_comps=30)


### Plotting histogram
scrub.plot_histogram()
plt.savefig(os.path.join(outdir,f"{os.path.basename(os.path.normpath(counts_matrix_dir))}_doublet_score_histogram.png"))

### Export
results = pd.Series(scrub.predicted_doublets_, name="scrublet_DropletType")
scores = pd.Series(scrub.doublet_scores_obs_, name="scrublet_Scores")
dataframe = pd.concat([barcodes_filtered_df, results, scores], axis=1)
dataframe.scrublet_DropletType = dataframe.scrublet_DropletType.replace(True, "doublet")
dataframe.scrublet_DropletType = dataframe.scrublet_DropletType.replace(False, "singlet")


print('Exporting results to {}'.format(os.path.join(outdir,f"{os.path.basename(os.path.normpath(counts_matrix_dir))}_doublet_calling_scrublet.csv")))


dataframe.to_csv( os.path.join(outdir,f"{os.path.basename(os.path.normpath(counts_matrix_dir))}_doublet_calling_scrublet.csv"), sep = ",", index = False)

