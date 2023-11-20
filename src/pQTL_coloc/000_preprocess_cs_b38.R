#!/usr/bin/env Rscript

#Rationale: Create separate credible sets file for each genomic locus.

args     = commandArgs(trailingOnly = TRUE)
cred_set = args[1]
tmp_path = args[2]
cred_set_b38 = args[3]

library(tidyverse)
library(data.table)

cs <- fread(cred_set)
b38 <- fread(cred_set_b38)
colnames(b38) <- c("chr_b38","pos_b38","pos1_b38")
cs_b38 <- cbind(cs,b38)

pheno <- unique(cs_b38$locus)

#split credset df into separate file for each locus:
##NB: allele2 is the effect allele, I change them because in /data/gen1/UKBiobank_500K/severe_asthma/Noemi_PhD/data/maf001_broad_pheno_1_5_ratio_betase_input_mungestat a1 is the effect allele, a2 is the non effect allele

for (i in pheno){
cs_b38_tmp <- cs_b38 %>% filter(locus == as.character(i))
locus <- as.character(unique(cs_b38_tmp$locus))
fwrite(cs_b38_tmp, paste0(tmp_path,"SA_",locus,"_b38"), quote=F, sep="\t")}