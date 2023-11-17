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
N = as.integer(mean(eqtlGWAS$N))
eqtlGWAS$N = N
eqtlGWAS <- mutate(eqtlGWAS,MAF=ifelse(freq>0.5,1-freq,freq))
eqtlGWAS <- mutate(eqtlGWAS,beta=ifelse(freq>0.5,-1*b,b))
eqtlGWAS <- eqtlGWAS %>% mutate(varbeta = se^2)
eqtlGWAS$allele1.eqtl <- toupper(eqtlGWAS$Allele1)
eqtlGWAS <- select(eqtlGWAS,"snp","MAF","beta","se","varbeta","pval","allele1.eqtl","N")

GWAS <- GWAS %>% filter(snp %in% eqtlGWAS$snp)

############################
# * if also want coloc results
############################
coloc_all = coloc.abf(dataset1=list(beta=GWAS$beta, varbeta=GWAS$varbeta,
                                    N=GWAS$N, type="cc", MAF=GWAS$MAF, snp=GWAS$snp),
                      dataset2=list(beta=eqtlGWAS$beta, varbeta=eqtlGWAS$varbeta,
                                    N=eqtlGWAS$N, type="quant", MAF=eqtlGWAS$MAF, snp=eqtlGWAS$snp))

saveRDS(coloc_all, paste0("/scratch/gen1/nnp5/Var_to_Gen_tmp/results/ubclung/", cred_set, "_ubclung_", probe, "_all_coloc.rds"))


########################
############################
# 1) Format LD matrix
############################
ld = data.table::fread(
  paste0(tmp_path, signal, "_", pheno, "_", tissue, ".raw")
) %>%
  select(!c(FID:PHENOTYPE))

ld %>%
  names %>%
  as_tibble %>%
  separate(value,c("c", "p", "a1", "a2")) %>%
  mutate(snp = paste0(c, "_", p, "_", pmin(a1, a2), "_", pmax(a1, a2))) %>%
  select(snp) %>%
  pull -> snps

ld = t(ld) %>%
  as_tibble %>%
  mutate(var=apply(., 1, var, na.rm=TRUE),
         snp=snps) %>%
  filter(var!=0) %>%
  filter(snp %in% GWAS$snp) %>%
  filter(snp %in% eqtlGWAS$snp)

N = ncol(ld)-2
X = t(as.matrix(ld[,1:N]))
colnames(X) = ld %>%
  select(snp) %>%
  pull # Sdd names to facilitate filtering Nas (still some NAs even with pairwise.complete.obs)
LDmatrix = cor(X, use="pairwise.complete.obs") # Pearson's correlation, takes some time...
LDmatrix[is.na(LDmatrix)] <- 0

############################
# 2) Format asthma GWAS
############################
bim = read_tsv(
    paste0(tmp_path, signal, "_", pheno, "_", tissue, ".bim"),
    col_names=c("chr", "snp", "morgans", "position", "allele1", "allele2")
  ) %>%
  select(!morgans) # Read bim file to check same set of SNPs and allele alignment

GWAS_df <- GWAS %>%
  filter(snp %in% colnames(LDmatrix)) %>% # Must be the same set of SNPs
  inner_join(bim) %>%
  mutate(beta=ifelse(allele1!=allele1.gwas, -(beta), beta)) # Check allele alignment

#Check GWAS has at least one variant with pval <= 0.000005:
GWAS_sign <- GWAS_df %>% filter(pval <= 0.000005)
if (dim(GWAS_sign)[1] < 1){
    stop(paste0(signal," ",tissue, " ", probe, " UBClung: No coloc.susie is possible because GWAS data has pvalue >= 5x10-6."))
}

GWAS_df %>%
  as.list -> GWAS

GWAS$type = "cc"
N = as.integer(mean(GWAS$N))
GWAS$N = N
GWAS$LD = LDmatrix

############################
# 3) Format eQTL GWAS
############################
eqtlGWAS_df <- eqtlGWAS  %>%
  filter(snp %in% colnames(LDmatrix)) %>% # Must be the same set of SNPs
  inner_join(bim) %>%
  mutate(beta=ifelse(allele1!=allele1.eqtl, -(beta), beta)) # Check allele alignment

#Check eqtlGWAS has at least one variant with pval <= 0.000005:
eqtlGWAS_sign <- eqtlGWAS_df %>% filter(pval <= 0.000005)
if (dim(eqtlGWAS_sign)[1] < 1){
    stop(paste0(signal," ",tissue, " ", probe, " UBClung: No coloc.susie is possible because eQTL-GWAS data has pvalue >= 5x10-6."))
} else {
    paste0(signal," ",tissue, " ", probe, " UBClung: Starting coloc.susie because both GWAS and eQTL data has pvalue <= 5x10-6.")
}

eqtlGWAS <- eqtlGWAS_df %>% as.list

eqtlGWAS$type = "quant"
eqtlGWAS$LD = LDmatrix
N = as.integer(mean(eqtlGWAS$N))
eqtlGWAS$N = N

############################
# Check datasets are ok
############################
check_dataset(eqtlGWAS, req="LD") -> eqtlGWAS_check
check_dataset(GWAS, req="LD") -> GWAS_check

if (is.null(GWAS_check) & is.null(eqtlGWAS_check)){
  cat(paste0(cred_set, "_", tissue, "_", probe, ": ok", "\n"),
      file=paste0(tmp_path,"GWAS_coloc_susie_check.txt"), append=TRUE)
} else {
  cat(paste0(cred_set, "_", tissue, "_", probe, ": warning", "\n"),
      file=paste0(tmp_path,"GWAS_coloc_susie_check.txt"), append=TRUE)
}

############################
## Run SuSie
############################
susie_GWAS = runsusie(GWAS, r2.prune=0.2, check_R=FALSE)
susie_eQTLGWAS = runsusie(eqtlGWAS, r2.prune=0.2, check_R=FALSE, verbose=TRUE)
susie_all = coloc.susie(susie_GWAS, susie_eQTLGWAS)

saveRDS(susie_all, paste0("/scratch/gen1/nnp5/Var_to_Gen_tmp/results/ubclung/", cred_set, "_", tissue, "_", probe, "_all_susie.rds"))
write_tsv(susie_all$summary %>% as_tibble,
          paste0("/scratch/gen1/nnp5/Var_to_Gen_tmp/results/ubclung/", cred_set, "_", tissue, "_", probe, "_susie_summary.tsv"))