#!/usr/bin/env Rscript

#Rationale: combine names in a document to create LD .raw file

#intermediate files in:
tmp_path="/scratch/gen1/nnp5/Var_to_Gen_tmp/"

library(tidyverse)


df <- data.frame(cred_set = c("SA_8_81292599_A_C","SA_6_90963614_A_AT","SA_5_110401872_T_C","SA_2_242692858_C_T",
        "SA_15_67442596_C_T","SA_12_56435504_G_C","SA_11_76296671_A_G","SA_9_6209697_A_G","SA_6_32586794_T_G",
        "SA_5_131885240_G_C","SA_3_33042712_T_C","SA_2_102913642_A_AAAAC","SA_17_38168828_A_G",
        "SA_16_27359021_A_C","SA_15_61068954_T_C","SA_12_48202941_C_T","SA_10_9064716_T_C")) %>%
  as_tibble %>%
  separate(cred_set, c("pheno", "chr", "pos", "a1", "a2"), sep="_")

df$ancestry <- "EUR"

df <- df %>%
  unite("snp", chr:a2, remove=FALSE) %>%
  mutate(eQTL  = "gtexlung",
         start = as.numeric(pos) - 5e5,
         end   = as.numeric(pos) + 5e5) %>%
  select(pheno, ancestry, snp, chr, pos, eQTL, start, end)

write_tsv(x=df, file=paste0(tmp_path,"gtex_Pairs_lookup.txt"),col_names=F)
