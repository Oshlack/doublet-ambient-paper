#!/usr/env/bin python

################################################################################
# RUN SCRIPT FOR DOUBLETDETECTION DOUBLET CALLING
################################################################################
# reference: https://github.com/jonathanshor/doubletdetection
# this script is adapted from demuxafy: https://demultiplexing-doublet-detecting-docs.readthedocs.io/en/latest/DoubletDetection.html

import os
import numpy as np
import doubletdetection
import tarfile
import matplotlib
matplotlib.use('PDF')
import matplotlib.pyplot as plt
import sys
import pandas as pd
import argparse

# Get path of mods directory from current script directory
mods_path = "/Users/meiermichelle/Library/CloudStorage/OneDrive-PeterMac/OSHLACK_LAB/PhD/ambient_doublet_singlecell/code/doublet_tools/" ## this is the path to the custom script for loading 10x data in python
sys.path.append(mods_path)
import read10x

## Set up parameters from trailing inputs
parser = argparse.ArgumentParser()

parser.add_argument("outdir", help="Output directory")
parser.add_argument("counts_matrix", help="Path to counts matrix directory")


args = parser.parse_args()

counts_matrix = args.counts_matrix
outdir = args.outdir


# need to make this function for multispawn reasons
def main():

    if not os.path.isdir(outdir):
        os.mkdir(outdir)

    ### Read in data ###
    raw_counts = read10x.import_cellranger_mtx(counts_matrix)

    try:
        barcodes_df = read10x.read_barcodes(counts_matrix + "/barcodes.tsv.gz")
    except:
        try:
            barcodes_df = read10x.read_barcodes(counts_matrix + "/barcodes.tsv")
        except:
            print("No barcode file found.")

    print('Counts matrix shape: {} rows, {} columns'.format(raw_counts.shape[0], raw_counts.shape[1]))

    # Remove columns with all 0s
    zero_genes = (np.sum(raw_counts, axis=0) == 0).A.ravel()
    raw_counts = raw_counts[:, ~zero_genes]

    print('Counts matrix shape after filtering: {} rows, {} columns'.format(raw_counts.shape[0], raw_counts.shape[1]))

    clf = doubletdetection.BoostClassifier( #https://doubletdetection.readthedocs.io/en/latest/tutorial.html
      n_iters=10, 
      clustering_algorithm="leiden", 
      clustering_kwargs =  {"flavor": "igraph"},
      standard_scaling=True,
      pseudocount=0.1,
      n_jobs=-1,
    )

    doublets = clf.fit(raw_counts).predict(p_thresh=1e-16, voter_thresh=0.5)

    results = pd.Series(doublets, name="DoubletDetection_DropletType")
    dataframe = pd.concat([barcodes_df, results], axis=1)

    dataframe.DoubletDetection_DropletType = dataframe.DoubletDetection_DropletType.replace(1.0, "doublet")
    dataframe.DoubletDetection_DropletType = dataframe.DoubletDetection_DropletType.replace(0.0, "singlet")

    outfile = os.path.join(
        outdir,
        f"{os.path.basename(os.path.normpath(counts_matrix))}_doublet_calling_DoubletDetection.csv"
    )

    print(f'Exporting results to {outfile}')
    dataframe.to_csv(outfile, sep=",", index=False)

    ### Figures ###
    print('Print threshold and convergence plots')

    doubletdetection.plot.convergence(
        clf,
        save=os.path.join(outdir, f"{os.path.basename(os.path.normpath(counts_matrix))}_convergence.pdf"),
        show=False,
        p_thresh=1e-16,
        voter_thresh=0.5
    )

    doubletdetection.plot.threshold(
        clf,
        save=os.path.join(outdir, f"{os.path.basename(os.path.normpath(counts_matrix))}_threshold.pdf"),
        show=False,
        p_step=6
    )


if __name__ == "__main__":
    main()
