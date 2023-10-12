library(data.table)
library(tidyverse)

args     = commandArgs(trailingOnly = TRUE)
cred_set = args[1]

tmp      = unlist(strsplit(cred_set, split="_"))
chr      = as.numeric(tmp[2])
location = as.numeric(tmp[3])

tmp_path="/scratch/gen1/nnp5/Var_to_Gen_tmp/"

gwas <- fread("/data/gen1/UKBiobank_500K/severe_asthma/Noemi_PhD/data/maf001_broad_pheno_1_5_ratio_betase_input_mungestat") %>%
  filter(b37chr==chr, between(bp, location-5e5, location+5e5))

write_delim(x=gwas, file=paste0(tmp_path, cred_set, "_GWASpairs.txt.gz"))