#!/usr/bin/env Rscript

#Rationale: merge genes from eQTL colocalisation analysis

library(tidyverse)
library(data.table)
library(readxl)

genes_raw <- "input/var2genes_raw.xlsx"
gtex <- read_excel(genes_raw, sheet = "GTExV8_eQTL_genes_symbol",col_names = c("ensembl","gene"))
eqtlgen <- read_excel(genes_raw, sheet = "eqtlGen_eQTL_genes",col_names = "gene")
ubclung <- read_excel(genes_raw, sheet = "UBClung_eQTL_genes",col_names = "gene")

gtex$ensembl <- NULL
eqtl_gene <- rbind(gtex,eqtlgen,ubclung) %>% unique()
eqtl_gene$evidence <- as.factor("eQTL")

write_tsv(x=eqtl_gene, file="input/eqtl_gene_merge")