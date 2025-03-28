#!/usr/bin/env Rscript

#Rationale: Collate results together and extract significant Variant/Locus-Tissue-Gene colocalisation, from UBCLung

library(plyr)
library(tidyverse)
library(purrr)
library(data.table)


#find probset that colocalised and map to gene using /data/gen1/reference/lung_eQTL/tabMerged_anno.txt

###############
### UBCLung ###
###############
print("UBCLung eQTL coloc results")
# .rds files containing standard coloc results with eQTLs
setwd("/scratch/gen1/nnp5/Var_to_Gen_tmp/results/ubclung/")
file_paths_ubclung = list.files(pattern="_all_coloc.rds")

# Manipulate file names to get various bits of info
file_paths_ubclung %>%
  as_tibble %>%
  mutate(file = gsub(x=value, pattern="_all_coloc.rds", replacement="")) %>%
  separate(col=file,
           into=c("pheno", "chr", "pos", "a1", "a2", "t1"),
           sep="_",
           remove=FALSE) -> file_split

# Clunky way to get out gene and tissue names from file names
tmp = strsplit(file_split$file, "_")
probe_1 = sapply(tmp, function(innerList) innerList[[length(innerList) - 1]])
probe_2 = sapply(tmp, tail, 1)
tissue_fct = function(df) {
  x = paste0(df[6])
}
tissue = sapply(tmp, tissue_fct)

# Add gene and tissue to above dataframe
file_split %>%
  mutate(probe_1   = probe_1,
         probe_2   = probe_2,
         tissue = tissue) %>%
  unite("snp", chr:a2) %>%
  select(file=value, pheno, snp, probe_1, probe_2, tissue) -> file_split

file_split$ProbeSet <- paste0(file_split$probe_1,"_",file_split$probe_2)

# Read all .rds files
data_list = lapply(file_paths_ubclung, readRDS)

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

#map probset to gene:
gene_probe <- fread("/data/gen1/reference/lung_eQTL/tabMerged_anno.txt",header=T)
gene_probe$ProbeSet <- gsub("_at","",gene_probe$ProbeSet)
setnames(gene_probe,"gene","gene_id")
data_summary <- left_join(data_summary, gene_probe,by="ProbeSet")

# Check if there are any colocalisations
data_summary %>% filter(coloc) %>% select(snp, pheno, gene_id)

# Save all results to file
#write_tsv(x=data_summary, file="/scratch/gen1/nnp5/Var_to_Gen_tmp/results/coloc_asthma_ubclung.tsv")
#for additional colocalisation analysis March 2025:
write_tsv(x=data_summary, file="/scratch/gen1/nnp5/Var_to_Gen_tmp/results/coloc_asthma_ubclung_addcredset_March2025.tsv")

# .rds files containing coloc.susie results with eQTLs
file_paths_ubclung = list.files(pattern="_all_susie.rds")

# Manipulate file names to get various bits of info
file_paths_ubclung %>%
  as_tibble %>%
  mutate(file = gsub(x=value, pattern="_all_susie.rds", replacement="")) %>%
  separate(col=file,
           into=c("pheno", "chr", "pos", "a1", "a2", "t1"),
           sep="_",
           remove=FALSE) -> file_split

# Clunky way to get out gene and tissue names from file names
tmp = strsplit(file_split$file, "_")
probe_1 = sapply(tmp, function(innerList) innerList[[length(innerList) - 1]])
probe_2 = sapply(tmp, tail, 1)
tissue_fct = function(df) {
  x = paste0(df[6])
}
tissue = sapply(tmp, tissue_fct)

# Add gene and tissue to above dataframe
file_split %>%
  mutate(probe_1   = probe_1,
         probe_2   = probe_2,
         tissue = tissue) %>%
  unite("snp", chr:a2) %>%
  select(file=value, pheno, snp, probe_1, probe_2, tissue) -> file_split
file_split$ProbeSet <- paste0(file_split$probe_1,"_",file_split$probe_2)

# Read all .rds files
data_list = lapply(file_paths_ubclung, readRDS)

# Function to extract coloc.susie information from results
summary_fct = function(df) {
  return(df$summary)
}

# Collapse list of coloc results into single dataframe and add a few new columns
##Additional commands to work around summary with NULL data:
tmp_1 <- lapply(data_list, summary_fct)
nonnull <- sapply(tmp_1, typeof)!="NULL"  # find all NULLs to omit
ID <- file_split$snp
index_row <- 1:length(tmp_1) #the length of tmp_1 which is the total nuber of all_susie.rds files
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

#map probset to gene:
data_summary_colocsusie <- left_join(data_summary_colocsusie, gene_probe,by="ProbeSet")

# Check if there are any colocalisations
data_summary_colocsusie %>% filter(coloc_susie) %>% select(snp, pheno, gene_id)

# Save all results to file
#write_tsv(x=data_summary_colocsusie, file="/scratch/gen1/nnp5/Var_to_Gen_tmp/results/colocsusie_asthma_ubclung.tsv")
#for additional colocalisation analysis March 2025:
write_tsv(x=data_summary_colocsusie, file="/scratch/gen1/nnp5/Var_to_Gen_tmp/results/colocsusie_asthma_ubclung_addcredset_March2025.tsv")

#SAVE COLOC AND COLOC.SUSIE GENES INTO XLSX FILE:
#coloc.susie does not have any significant colocalisation:
gene_coloc <- as.data.frame(data_summary %>% filter(coloc) %>% select(gene_id))
gene_coloc_susie <- as.data.frame(data_summary_colocsusie %>% filter(coloc_susie) %>% select(gene_id))
print(gene_coloc)
#NO gene coloc_susie
write.table(gene_coloc,"/home/n/nnp5/PhD/PhD_project/Var_to_Gene/input/UBCLungeqtl_var2genes_raw.txt", row.names=FALSE, col.names=FALSE, quote=F)



#then add the genes in the /Var_to_Gene/input/var2genes_raw.xlsx file
