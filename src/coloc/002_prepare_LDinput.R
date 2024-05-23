#!/usr/bin/env Rscript

#Rationale: combine names in a document to create LD .raw file

args     = commandArgs(trailingOnly = TRUE)
tissue = args[1]

#intermediate files in:
tmp_path="/scratch/gen1/nnp5/Var_to_Gen_tmp/"

library(tidyverse)


#df <- data.frame(cred_set = c("SA_8_81292599_C_A","SA_6_90963614_AT_A","SA_5_110401872_C_T",
#    "SA_2_242692858_T_C","SA_15_67442596_T_C","SA_12_56435504_C_G","SA_11_76296671_G_A",
#    "SA_9_6209697_G_A","SA_5_131885240_C_G","SA_3_33042712_C_T",
#    "SA_2_102913642_AAAAC_A","SA_17_38168828_G_A","SA_16_27359021_C_A","SA_15_61068954_C_T",
#    "SA_12_48202941_T_C","SA_10_9064716_C_T")) %>%
#  as_tibble %>%
#  separate(cred_set, c("pheno", "chr", "pos", "a1", "a2"), sep="_")

df <- data.frame(cred_set = "SA_3_50024027_C_CA") %>%
  as_tibble %>%
  separate(cred_set, c("pheno", "chr", "pos", "a1", "a2"), sep="_")

df$ancestry <- "EUR"

df <- df %>%
  unite("snp", chr:a2, remove=FALSE) %>%
  mutate(eQTL  = tissue,
         start = as.numeric(pos) - 5e5,
         end   = as.numeric(pos) + 5e5) %>%
  select(pheno, ancestry, snp, chr, pos, eQTL, start, end)

write_tsv(x=df, file=paste0(tmp_path,"gtex_Pairs_lookup.txt"),col_names=F, append=T)
