#!/usr/bin/env Rscript

################################################################################
# RUN SCRIPT FOR SCDBLFINDER DOUBLET CALLING - WITH FILTERED BARCODES
################################################################################
# reference: https://github.com/plger/scDblFinder


### load libraries ----
suppressPackageStartupMessages({
  library(scDblFinder)
  library(SingleCellExperiment)
  library(tidyverse)
  library(readr)
  library(here)
  library(Seurat)
  library(scater)
})

### input parameters ----
args <- commandArgs(trailingOnly = TRUE)

# Check that we get two input arguments in right order
if (length(args) != 2) {
  stop("Expecting two input arguments: output directory, input counts", call. = FALSE)
}

# take from args
output_dir <- args[1]
counts_dir <- args[2]

message(sprintf("arguments detected:\noutput:%s\ncounts:%s", output_dir, counts_dir))

### Read data ----
counts <- Read10X(as.character(counts_dir), gene.column = 1)
sce <- SingleCellExperiment(assays = c("counts" = counts))
#filter based on barcodes
barcodes <- read_tsv(paste0(counts_dir, "/filtered_qc_barcodes.tsv"), col_names = T, show_col_types = F) %>% pull(".")
sce <- sce[, colnames(sce) %in% barcodes]

### Doublet calling ----

#run pre-processing required
message("Run PCA and UMAP")
set.seed(1997)
sce <- logNormCounts(sce) %>%
  runPCA() %>%
  runUMAP()


message("Run scDblFinder")
set.seed(1997)
sce <- scDblFinder(sce,
                   clusters = T,
                   dbr.sd = 1) #no assumption about doublet rates, as recommended by scDblFinder people

gc()



### Make dataframe ----
#inspired by: https://demultiplexing-doublet-detecting-docs.readthedocs.io/en/latest/Scds.html
Doublets <- as.data.frame(cbind(rownames(colData(sce)), colData(sce)$scDblFinder.score, as.character(colData(sce)$scDblFinder.class)))
colnames(Doublets) <- c("Barcode","scDblFinder_score","scDblFinder_DropletType")

### Export ----
filename <- sprintf("%s/%s_doublet_calling_scDblFinder.csv", output_dir, basename(counts_dir))
message(sprintf("Writing results to %s", filename))
write_csv(x = Doublets, file = filename)

sessionInfo()

