#!/usr/bin/env Rscript

################################################################################
# RUN SCRIPT FOR DOUBLETFINDER DOUBLET CALLING
################################################################################
# reference: https://github.com/chris-mcginnis-ucsf/DoubletFinder


### load libraries ----
suppressPackageStartupMessages({
  library(DoubletFinder)
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
if (length(args) != 3) {
  stop("Expecting three input arguments: output directory, input counts, expected #doublets", call. = FALSE)
}

# take from args
output_dir <- args[1]
counts_dir <- args[2]
n_doublets <- as.numeric(args[3])

message(sprintf("arguments detected:\noutput:%s\ncounts:%s\nndoublets:%i", output_dir, counts_dir, n_doublets))

### Read data ----
counts <- Read10X(as.character(counts_dir), gene.column = 1)
seu <- CreateSeuratObject(counts = counts)

### sc-wflow ----

#doubletfinder requires normalisation etc.
seu <- NormalizeData(seu, verbose = F) %>%
  FindVariableFeatures(., verbose = F) %>%
  ScaleData(., verbose = F) %>%
  RunPCA(., verbose = F) %>%
  RunUMAP(., dims = 1:30, verbose = F)


### Parameter sweep ----
# inspired by https://github.com/drneavin/Demultiplexing_Doublet_Detecting_Docs/blob/main/scripts/DoubletFinder.R
sweep.res.list <- paramSweep(seu, PCs = 1:30, sct = FALSE) #no SCT was performed
sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE) #GT= ground truth, not being used here
bcmvn <- find.pK(sweep.stats)
pK_select <- bcmvn$pK[which(bcmvn$BCmetric == max(bcmvn$BCmetric))] %>% as.character() %>% as.numeric() #automate this, get highest one

#save plot
plot <- ggplot(bcmvn, aes(pK, BCmetric, group=1)) +
  geom_point() +
  geom_line() +
  geom_vline(xintercept = as.factor(pK_select), linetype="dashed",
             color = "red", size=1)+
  theme_bw()
ggsave(plot, filename = sprintf("%s/pKvBCmetric.png", output_dir), width = 10, height = 6)


### Doublet calling ----
# using input n_doublets -> not adjusting for hom typic doublets (as they describe in their vignette)
seu <- doubletFinder(seu, PCs = 1:30, pN = 0.25, pK =pK_select, nExp = n_doublets, reuse.pANN = NULL, sct = FALSE) #pN= default

gc()


### Make dataframe ----
doublets <- as.data.frame(cbind(colnames(seu),
                                seu@meta.data[,paste0("pANN_0.25_",pK_select, "_", n_doublets)],
                                seu@meta.data[,paste0("DF.classifications_0.25_",pK_select, "_", n_doublets)]))

colnames(doublets) <-  c("Barcode","DoubletFinder_score","DoubletFinder_DropletType")
doublets$DoubletFinder_DropletType <- gsub("Singlet","singlet",doublets$DoubletFinder_DropletType) %>% gsub("Doublet","doublet",.)



### Export ----
filename <- sprintf("%s/%s_doublet_calling_DoubletFinder.csv", output_dir, basename(counts_dir))
message(sprintf("Writing results to %s", filename))
write_csv(x = doublets, file = filename)

sessionInfo()

