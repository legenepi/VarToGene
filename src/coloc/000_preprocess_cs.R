#!/usr/bin/env Rscript

#Rationale: Create separate credible sets file for each genomic locus.

args     = commandArgs(trailingOnly = TRUE)
cred_set = args[1]
tmp_path = args[2]

library(tidyverse)
library(data.table)

cs <- fread(cred_set)
#filter out the locus on the MHC region
pheno_tmp <- cs %>% filter(locus != "6_rs9271365_32086794_33086794")
pheno <- unique(pheno_tmp$locus)

#split credset df into separate file for each locus:
##NB: allele2 is the effect allele, I change them because in /data/gen1/UKBiobank_500K/severe_asthma/Noemi_PhD/data/maf001_broad_pheno_1_5_ratio_betase_input_mungestat a1 is the effect allele, a2 is the non effect allele

for (i in pheno){
cs_tmp <- cs %>% filter(locus == as.character(i))
chr <- cs_tmp %>% filter(PIP_average == max(cs_tmp$PIP_average))  %>% select(chromosome)
location <- cs_tmp %>% filter(PIP_average == max(cs_tmp$PIP_average)) %>% select(position)
allele1 <- cs_tmp %>% filter(PIP_average == max(cs_tmp$PIP_average))  %>% select(allele2)
allele2 <- cs_tmp %>% filter(PIP_average == max(cs_tmp$PIP_average))  %>% select(allele1)
fwrite(cs_tmp, paste0(tmp_path,"SA_",chr,"_",location,"_",allele1, "_",allele2), quote=F, sep="\t")}

