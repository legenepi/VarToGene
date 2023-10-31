#!/usr/bin/env Rscript

#==================================================================
# File:        000A_eqtl_gtex_extraction.R
# Project:     SevereAsthma
# Author:      NNP - edited from KC
# Date:        31 October 2023
# Rationale:   Load eQTL data for relevant tissues and and extract
#              chromosome and positions for liftOver.
#==================================================================

suppressMessages(library(tidyverse))
suppressMessages(library(data.table))
#Sys.setenv(LIBARROW_BINARY = TRUE); install.packages('arrow', type = "source")
suppressMessages(library(arrow))

options(scipen=999)

args <- commandArgs(T)

chr <- args[1]

setwd("/data/gen1/reference/GTEx_Analysis_v8_QTLs/GTEx_Analysis_v8_EUR_eQTL_all_associations")

# Identify files ----
file.names = list.files(path = ".", pattern = paste0("*.v8.EUR.allpairs.chr", chr, ".parquet")) %>%
  #Kayesha:
  #str_subset(pattern = "Adrenal|Artery|Blood|Brain|Esophagus|Heart|Ileum|Liver|Lung|Muscle|Nerve|Pituitary|Spleen|Stomach")
  str_subset(pattern = "Skin|Colon")

for (file in file.names) {

  # Read in data ----
  name = str_remove(file, ".parquet")
  cat(paste0("Tissue: ", name, "\n"))
 
  cat("Loading data...\n")
  df <- read_parquet(file)

  # Separate variant_id column into chromosome, position, reference allele and alternative allele ----
  cat("Wrangling data...\n")
  sep <- df %>%
    separate(variant_id, c("CHROM", "POS", "REF", "ALT")) %>%
    mutate(CHROM = str_replace(CHROM, "chr", ""),
           CHROM = str_replace(CHROM, "X", "23"),
           CHROM = as.numeric(CHROM),
           POS = as.numeric(POS),
           ID = paste0(CHROM, "_", POS, "_", pmin(REF, ALT), "_", pmax(REF, ALT)))

    ## Write file ----
    cat(paste0("Writing file: ", name, ".hg38.txt.gz", "\n"))
    fwrite(sep, file = paste0("/scratch/gen1/nnp5/Var_to_Gen_tmp/liftover_gtexv8/", name, ".hg38.txt.gz"), quote = F, sep = "\t", na = "NA", row.names = F, col.names = T, compress = "auto")


  # Generate BED file for liftOver ----
  cat(paste0("Creating bed file: ", name, ".hg38.bed", "\n"))
  sep %>%
    select(chrom = CHROM, chromStart = POS, ID) %>%
    mutate(chromEnd = as.integer(chromStart+1),
           chrom = paste0("chr", chrom),
           chrom = str_replace(chrom, "chr23", "chrX")) %>%
    select(chrom, chromStart, chromEnd, ID) %>%
    fwrite(file = paste0("/scratch/gen1/nnp5/Var_to_Gen_tmp/liftover_gtexv8/bed/", name, ".hg38.bed"), quote = F, sep = "\t", na = "NA", row.names = F, col.names = F)
  
  cat("Done.\n\n")
  
  }
