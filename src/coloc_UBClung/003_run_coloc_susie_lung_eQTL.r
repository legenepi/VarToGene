#!/usr/bin/env Rscript

#Rationale: Run coloc and coloc.susie in UBCLung eQTL

suppressMessages(library(tidyverse))
library(coloc)
suppressMessages(library(Rfast))
library(data.table)

args     = commandArgs(trailingOnly = TRUE)
tissue   = args[1]
cred_set = args[2]
probe    = args[3]

tmp_path="/scratch/gen1/nnp5/Var_to_Gen_tmp/ubclung/"

cred_set %>%
  as_tibble %>%
  separate(value, c("pheno", "chr", "pos", "a1", "a2"), sep="_") %>%
  select(pheno) %>%
  pull -> pheno

cred_set %>%
  as_tibble %>%
  separate(value, c("pheno", "chr", "pos", "a1", "a2"), sep="_") %>%
  unite("signal", chr:a2) %>%
  pull -> signal

##allele are place in alphabetical order ! (no dependency on effect allele, reference/alternative)

############################
## Read asthma GWAS sumstat, eQTL GWAS, and LD matrix
############################
GWAS = read_delim(paste0("/scratch/gen1/nnp5/Var_to_Gen_tmp/", cred_set, "_GWASpairs.txt.gz")) %>%
  rename(chr = b37chr, pos = bp, beta = LOG_ODDS, SE = se) %>%
  select(-snpid) %>%
  arrange(chr, pos) %>%
  drop_na(beta) %>%
  mutate(varbeta = SE^2) %>%
  rename(position=pos)
GWAS$allele1.gwas <- GWAS$a1
GWAS <- GWAS %>% mutate(snp = paste0(chr, "_", position, "_", pmin(a1, a2), "_", pmax(a1, a2))) %>%
  select(c(snp,beta,SE,eaf,pval,MAF,varbeta,allele1.gwas)) %>%
  distinct(snp, .keep_all=TRUE)

GWAS$N <- as.numeric(38405)

eQTL <- fread(paste0(tmp_path,"eQTL_region_stat/",cred_set,"_eQTL_region.txt",sep=""),header=T)
setnames(eQTL,"#Probe","ProbeSet")
eQTL <- eQTL %>% filter(ProbeSet==probe)
eQTL <- mutate(eQTL,snp=paste(CHR,BP,pmin(toupper(Allele1),toupper(Allele2)),pmax(toupper(Allele1),toupper(Allele2)),sep="_"))
eqtlGWAS <- eQTL %>% filter(snp %in% GWAS$snp)

setnames(eqtlGWAS,"Freq1","freq")
setnames(eqtlGWAS,"Effect","b")
setnames(eqtlGWAS,"StdErr","se")
setnames(eqtlGWAS,"P","pval")

eqtlGWAS$freq <- as.numeric(eqtlGWAS$freq)
eqtlGWAS$b <- as.numeric(eqtlGWAS$b)
eqtlGWAS$se <- as.numeric(eqtlGWAS$se)
eqtlGWAS$z <- eqtlGWAS$b/eqtlGWAS$se
eqtlGWAS$N <- 1/(2*eqtlGWAS$freq*(1-eqtlGWAS$freq)*eqtlGWAS$se*eqtlGWAS$se)-eqtlGWAS$z*eqtlGWAS$z
eqtlGWAS <- mutate(eqtlGWAS,MAF=ifelse(freq>0.5,1-freq,freq))
eqtlGWAS <- mutate(eqtlGWAS,beta=ifelse(freq>0.5,-1*b,b))
eqtlGWAS <- eqtlGWAS %>% mutate(varbeta = se^2)
eqtlGWAS <- select(eqtlGWAS,"snp","MAF","beta","se","varbeta","pval","N")

GWAS <- GWAS %>% filter(snp %in% eQTL$snp)

############################
# * if also want coloc results
############################
coloc_all = coloc.abf(dataset1=list(beta=GWAS$beta, varbeta=GWAS$varbeta,
                                    N=GWAS$N, type="cc", MAF=GWAS$MAF, snp=GWAS$snp),
                      dataset2=list(beta=eqtlGWAS$beta, varbeta=eqtlGWAS$varbeta,
                                    N=eqtlGWAS$N, type="quant", MAF=eqtlGWAS$MAF, snp=eqtlGWAS$snp))

saveRDS(coloc_all, paste0("/scratch/gen1/nnp5/Var_to_Gen_tmp/results/ubclung/", cred_set, "_eqtlGenWB_", probe, "_all_coloc.rds"))



########################
############################
# 1) Format LD matrix
############################
#TO BE done !!!