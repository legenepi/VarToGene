#!/usr/bin/env Rscript

#Rationale: Collate results together and extract significant Variant/Locus-Tissue-Gene colocalisation, from UBCLung

library(plyr)
library(tidyverse)
library(purrr)
library(xlsx)

TO BE EDITED FOR UBCLUNG !

find probset that colocalised and map to gene using /data/gen1/reference/lung_eQTL/tabMerged_anno.txt

#gene_probe <- fread("/data/gen1/reference/lung_eQTL/tabMerged_anno.txt",header=T)
#gene_probe$ProbeSet <- gsub("_at","",gene_probe$ProbeSet)
#setnames(gene_probe,"gene","gene_id")

###############
### UBCLung ###
###############
print("eqtlGen eQTL coloc results")
# .rds files containing standard coloc results with eQTLs
setwd("/scratch/gen1/nnp5/Var_to_Gen_tmp/results/eqtlgen/")
file_paths_eqtlgen = list.files(pattern="_all_coloc.rds")

# Manipulate file names to get various bits of info
file_paths_eqtlgen %>%
  as_tibble %>%
  mutate(file = gsub(x=value, pattern="_all_coloc.rds", replacement="")) %>%
  separate(col=file,
           into=c("pheno", "chr", "pos", "a1", "a2", "t1"),
           sep="_",
           remove=FALSE) -> file_split

# Clunky way to get out gene and tissue names from file names
tmp = strsplit(file_split$file, "_")
gene = sapply(tmp, tail, 1)
tissue_fct = function(df) {
  x = paste0(df[6])
}
tissue = sapply(tmp, tissue_fct)

# Add gene and tissue to above dataframe
file_split %>%
  mutate(gene   = gene,
         tissue = tissue) %>%
  unite("snp", chr:a2) %>%
  select(file=value, pheno, snp, gene, tissue) -> file_split

# Read all .rds files
data_list = lapply(file_paths_eqtlgen, readRDS)

# Function to extract coloc information from results
summary_fct = function(df) {
  return(df$summary)
}

# Collapse list of coloc results into single dataframe and add a few new columns
data_summary = ldply(lapply(data_list, summary_fct)) %>%
  as_tibble %>%
  bind_cols(file_split) %>%
  select(-file) %>%
  mutate(coloc = if_else(PP.H4.abf >= 0.9, TRUE, FALSE)) %>%
  arrange(desc(PP.H4.abf))

# Check if there are any colocalisations
data_summary %>% filter(coloc) %>% select(snp, pheno, gene)

# Save all results to file
write_tsv(x=data_summary, file="/scratch/gen1/nnp5/Var_to_Gen_tmp/results/coloc_asthma_eqtlgen.tsv")

# .rds files containing coloc.susie results with eQTLs
file_paths_gtex = list.files(pattern="_all_susie.rds")

# Manipulate file names to get various bits of info
file_paths_gtex %>%
  as_tibble %>%
  mutate(file = gsub(x=value, pattern="_all_susie.rds", replacement="")) %>%
  separate(col=file,
           into=c("pheno", "chr", "pos", "a1", "a2", "t1"),
           sep="_",
           remove=FALSE) -> file_split

# Clunky way to get out gene and tissue names from file names
tmp = strsplit(file_split$file, "_")
gene = sapply(tmp, tail, 1)
tissue_fct = function(df) {
  x = paste0(df[6])
}
tissue = sapply(tmp, tissue_fct)

# Add gene and tissue to above dataframe
file_split %>%
  mutate(gene   = gene,
         tissue = tissue) %>%
  unite("snp", chr:a2) %>%
  select(file=value, pheno, snp, gene, tissue) -> file_split

# Read all .rds files
data_list = lapply(file_paths_gtex, readRDS)

# Function to extract coloc.susie information from results
summary_fct = function(df) {
  return(df$summary)
}

# Collapse list of coloc results into single dataframe and add a few new columns
##Additional commands to work around summary with NULL data:
tmp_1 <- lapply(data_list, summary_fct)
nonnull <- sapply(tmp_1, typeof)!="NULL"  # find all NULLs to omit
ID <- file_split$snp
index_row <- 1:169 #the length of tmp_1 which is the total nuber of all_susie.rds files
tmp_2 <- map2(tmp_1, ID, ~cbind(.x, snp = .y))
tmp_3 <- map2(tmp_2, index_row, ~cbind(.x, n_index = .y))
tmp_4 <- tmp_3[nonnull]

file_split$n_index <- index_row

data_summary_colocsusie = ldply(tmp_4) %>%
  as_tibble %>%
  left_join(file_split,by=c("snp","n_index")) %>%
  select(-file) %>%
  mutate(coloc_susie = if_else(PP.H4.abf >= 0.9, TRUE, FALSE)) %>%
  arrange(desc(PP.H4.abf))

# Check if there are any colocalisations
data_summary_colocsusie %>% filter(coloc_susie) %>% select(snp, pheno, gene)


# Save all results to file
write_tsv(x=data_summary_colocsusie, file="/scratch/gen1/nnp5/Var_to_Gen_tmp/results/colocsusie_asthma_eqtlgen.tsv")

#SAVE COLOC AND COLOC.SUSIE GENES INTO XLSX FILE:
#coloc.susie does not have any significant colocalisation:
gene_coloc <- as.data.frame(data_summary %>% filter(coloc) %>% select(gene))
write.xlsx(gene_coloc,"/home/n/nnp5/PhD/PhD_project/Var_to_Gene/input/var2genes_raw.xlsx", sheetName = "eqtlGen_eQTL_genes", row.names=FALSE, col.names=FALSE, append=TRUE)
