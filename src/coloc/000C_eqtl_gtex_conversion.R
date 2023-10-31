#!/usr/bin/env Rscript

#==================================================================
# File:        000C_eqtl_gtex_conversion.R
# Project:     SevereAsthma
# Author:      NNP-edited from KC
# Date:        31 October 2023
# Ratinale:    Update coordinates in original GTeX files to GRCh37.
#==================================================================

suppressMessages(library(tidyverse))
suppressMessages(library(data.table))

options(scipen=999)

args <- commandArgs(T)

chr <- args[1]

setwd("/scratch/gen1/nnp5/Var_to_Gen_tmp/liftover_gtexv8")

# Identify files ----
file.names = list.files(path = ".", pattern = paste0("*.v8.EUR.allpairs.chr", chr, ".hg38.txt.gz")) %>%
  #Kayesha
  #str_subset(pattern = "Adrenal|Artery|Blood|Brain|Esophagus|Heart|Ileum|Liver|Lung|Muscle|Nerve|Pituitary|Spleen|Stomach")
  str_subset(pattern = "Colon|Skin")

for (file in file.names) {


  # Read in data ----
  name = str_remove(file, ".hg38.txt.gz")
  cat(paste0("Tissue: ", name, "\n"))

  ## GTeX
  hg38 <- as_tibble(fread(file))
  ## Liftover file (hg19)
  liftover <- as_tibble(fread(paste0("bed/", name, ".hg19.bed"), header = FALSE))


  # Remove chromosomes from liftover bed (hg19) which are alternative contigs, and filter for only chromosome of GTeX data ----
  ## As multiple duplicate variants are present, remove these
  cat("Wrangling liftOver data...\n")
  if (chr == "X") {

    liftover <- liftover %>%
      select(CHROM = V1, GENPOS = V2, ID = V4) %>%
      mutate(CHROM = str_replace(CHROM, "chr", ""),
             CHROM = str_replace(CHROM, "X", "23")) %>%
      filter(CHROM == 23) %>%
      mutate(CHROM = as.integer(CHROM)) %>%
      distinct(ID, .keep_all = TRUE)
    #arrange(CHROM, GENPOS)

  } else {
  
    liftover <- liftover %>%
      select(CHROM = V1, GENPOS = V2, ID = V4) %>%
      mutate(CHROM = str_replace(CHROM, "chr", ""),
             CHROM = str_replace(CHROM, "X", "23")) %>%
      filter(CHROM == chr) %>%
      mutate(CHROM = as.integer(CHROM)) %>%
      distinct(ID, .keep_all = TRUE)
    }
  

  # Perform mapping ----
  cat("Performing mapping...\n")
  hg19 <- hg38 %>%
    right_join(liftover, by = "ID") %>%
    select(-CHROM.x, -POS) %>%
    rename(CHROM = CHROM.y, POS = GENPOS) %>%
    mutate(ID = paste0(CHROM, "_", POS, "_", pmin(REF, ALT), "_", pmax(REF, ALT))) %>%
    relocate(c(CHROM, POS, ID), .after = phenotype_id) %>%
    select(-tss_distance)        # Will no longer apply to new coordinates


  # Write file ----
  cat(paste0("Writing file: ", name, ".hg19.txt.gz", "\n"))
  fwrite(hg19, paste0(name, ".hg19.txt.gz"), quote = F, sep = "\t", na = "NA", row.names = F, col.names = T, compress = "auto")
  cat("Done.\n\n")
  }