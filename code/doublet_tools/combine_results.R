#!/usr/bin/env Rscript

################################################################################
# RUN SCRIPT TO COMBINE DOUBLET CALLING OUTCOMES
################################################################################

### load libraries ----
suppressPackageStartupMessages({
  library(tidyverse)
  library(readr)
})

### input parameters ----
args <- commandArgs(trailingOnly = TRUE)

# Check that we get two input arguments in right order
if (length(args) != 1) {
  stop("Expecting one input argument: output directory", call. = FALSE)
}

# take from args
output_dir <- args[1]

message(sprintf("arguments detected:\noutput:%s", output_dir))

### Read data ----
files <- list.files(path = output_dir, pattern = "doublet_calling.*\\.csv", full.names = T)
#remove consensus file, in case it exists and we want to overwrite it
files <- files[!grepl("consensus", files)]

all_results <- lapply(files, function(file){
  read_csv(file=file, show_col_types = FALSE)
  })
all_results <- lapply(all_results, function(x){
  x %>% column_to_rownames("Barcode") %>%   mutate(Barcode = rownames(.))
})

doublet_DF <- Reduce(full_join, all_results) %>%
  mutate(consensus_singlet = if_else(
           if_all(ends_with("DropletType"), ~ . == "singlet"), "singlet", "non-singlet"))

#reorder results based on barcodes file
barcodes <- read_tsv(paste0(dirname(output_dir), "/barcodes.tsv.gz"), show_col_types = F, col_names = F) %>% pull(X1)

doublet_DF_reordered <- doublet_DF[match(barcodes, doublet_DF$Barcode),]
#test
all(doublet_DF_reordered$Barcode == barcodes)

### Export ----
name <- gsub("_doublet_calling.*", "", basename(list.files(path = output_dir, pattern = "doublet_calling.*\\.csv", full.names = T)[1]))
filename <- sprintf("%s/%s_doublet_calling_consensus.csv", output_dir, name)
message(sprintf("Writing results to %s", filename))
write_csv(x = doublet_DF, file = filename)

sessionInfo()

