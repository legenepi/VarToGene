#!/usr/bin/env Rscript

#Rationale: Lookup of overlapping region between ssevere asthma GWAS and UBC Lung eQTL significant genes in order to
#perform, colocalisation only in regions that have significant associations in both studies.

sink(stderr())

suppressMessages(library(data.table))
suppressMessages(library(tidyverse))
suppressMessages(library(dplyr))

args     = commandArgs(trailingOnly = TRUE)
tissue   = args[1] #"UBCLung"
cred_set = args[2]

tmp      = unlist(strsplit(cred_set, split="_"))
chr_sentinel      = as.numeric(tmp[2])
location_sentinel = as.numeric(tmp[3])

tmp_path="/scratch/gen1/nnp5/Var_to_Gen_tmp/ubclung/"

##allele are place in alphabetical order ! (no dependency on effect allele, reference/alternative)

eGene <- fread("/data/gen1/LF_HRC_transethnic/eQTL/lung_eQTL/eGene_list.txt")
eGene$chr <- gsub("_[0-9]*_[A-Z]*_[A-Z]*","",eGene$SNP)
eGene$bp <- gsub("^[0-9]*_","",eGene$SNP)
eGene$bp <- gsub("_[A-Z]*_[A-Z]*","",eGene$bp)

eGENE_lookup <- eGene %>% filter(chr == chr_sentinel,
                                  bp >= location_sentinel-1000000,
                                  bp <= location_sentinel+1000000)
eGENE_lookup <-  eGENE_lookup %>% mutate(sentinel=cred_set)
probeset <- eGENE_lookup %>% distinct(ProbeSet)

#colnames: SNP ProbeSet p_value BF TESTS chr bp sentinel
write_delim(x=eGENE_lookup, file=paste0(tmp_path, cred_set,"_", tissue, "_pairs.txt.gz"))
write_delim(x=probeset, file=paste0(tmp_path, cred_set,"_", tissue, "_probesets.txt"), col_names=FALSE)
