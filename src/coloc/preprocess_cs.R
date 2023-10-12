#!/usr/bin/env Rscript

#Rationale: Create separate credible sets file for each genomic locus.

args     = commandArgs(trailingOnly = TRUE)
cred_set = args[1]
tmp_path = args[2]

library(tidyverse)
library(data.table)

cs <- fread(cred_set)
pheno <- unique(cs$locus)

#split credset df into separate file foe each locus:
for (i in pheno){
cs_tmp <- cs %>% filter(locus == as.character(i))
chr <- cs_tmp %>% filter(PIP_average == max(cs_tmp$PIP_average))  %>% select(chromosome)
location <- cs_tmp %>% filter(PIP_average == max(cs_tmp$PIP_average)) %>% select(position)
allele1 <- cs_tmp %>% filter(PIP_average == max(cs_tmp$PIP_average))  %>% select(allele1)
allele2 <- cs_tmp %>% filter(PIP_average == max(cs_tmp$PIP_average))  %>% select(allele2)
fwrite(cs_tmp, paste0(tmp_path,"SA_",chr,"_",location,"_",allele1, "_",allele2), quote=F, sep="\t")}

