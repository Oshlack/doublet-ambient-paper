#!/usr/bin/env Rscript

################################################################################
# RUN SCRIPT FOR SCDS DOUBLET CALLING
################################################################################
# reference: https://github.com/kostkalab/scds


### load libraries ----
suppressPackageStartupMessages({
  library(scds)
  library(SingleCellExperiment)
  library(tidyverse)
  library(readr)
  library(here)
  library(Seurat)
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


### Doublet calling ----

#Annotate doublet using co-expression based doublet scoring:
message("Start doublet calling: cxds")
sce <- cxds(sce, retRes=T, estNdbl=T)
gc()

#Annotate doublet using binary classification based doublet scoring:
message("Start doublet calling: bcds")
sce <- bcds(sce, retRes=T, estNdbl=T)
gc()

#Combine both annotations into a hybrid annotation
message("Start doublet calling: hybrid")
sce <- cxds_bcds_hybrid(sce, estNdbl=T)
gc()

### Make dataframe ----
#inspired by: https://demultiplexing-doublet-detecting-docs.readthedocs.io/en/latest/Scds.html
Doublets <- as.data.frame(cbind(rownames(colData(sce)), colData(sce)$hybrid_score, colData(sce)$hybrid_call))
colnames(Doublets) <- c("Barcode","scds_score","scds_DropletType")
Doublets$scds_DropletType <- gsub("FALSE","singlet",Doublets$scds_DropletType)
Doublets$scds_DropletType <- gsub("TRUE","doublet",Doublets$scds_DropletType)

### Export ----
filename <- sprintf("%s/%s_doublet_calling_scds.csv", output_dir, basename(counts_dir))
message(sprintf("Writing results to %s", filename))
write_csv(x = Doublets, file = filename)

sessionInfo()

