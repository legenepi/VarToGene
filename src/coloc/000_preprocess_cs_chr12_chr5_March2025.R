#!/usr/bin/env Rscript

#Rationale: Create separate credible sets file for each genomic locus.

args     = commandArgs(trailingOnly = TRUE)
cred_set = args[1]
tmp_path = args[2]

library(tidyverse)
library(data.table)

cs <- fread(cred_set)
cs <- cs %>% rename(chrom = Chr, position = Pos, snpid = SNP_ID)
locus <- unique(cs$Replicated_locus)

#split credset df into separate file for each locus:
##NB: allele2 is the effect allele, I change them because in /data/gen1/UKBiobank_500K/severe_asthma/Noemi_PhD/data/maf001_broad_pheno_1_5_ratio_betase_input_mungestat a1 is the effect allele, a2 is the non effect allele

for (i in locus){
cs_tmp <- cs %>% filter(Credible_set == 1)
cs_tmp <- cs_tmp %>% filter(Replicated_locus == as.character(i))
chrom <- cs_tmp %>% filter(PIP == max(cs_tmp$PIP))  %>% select(chrom)
location <- cs_tmp %>% filter(PIP == max(cs_tmp$PIP)) %>% select(position)
allele1 <- cs_tmp %>% filter(PIP == max(cs_tmp$PIP))  %>% select(allele2)
allele2 <- cs_tmp %>% filter(PIP == max(cs_tmp$PIP))  %>% select(allele1)
fwrite(cs_tmp, paste0(tmp_path,"SA_",chrom,"_",location,"_",allele1, "_",allele2), quote=F, sep="\t")}

for (i in locus){
cs_tmp <- cs %>% filter(Replicated_locus == as.character(i))
cs_tmp <- cs_tmp %>% filter(Credible_set == 2)
if (dim(cs_tmp)[1] > 0) {
chrom <- cs_tmp %>% filter(PIP == max(cs_tmp$PIP))  %>% select(chrom)
location <- cs_tmp %>% filter(PIP == max(cs_tmp$PIP)) %>% select(position)
allele1 <- cs_tmp %>% filter(PIP == max(cs_tmp$PIP))  %>% select(allele2)
allele2 <- cs_tmp %>% filter(PIP == max(cs_tmp$PIP))  %>% select(allele1)
fwrite(cs_tmp, paste0(tmp_path,"SA_",chrom,"_",location,"_",allele1, "_",allele2), quote=F, sep="\t")}}

warnings()