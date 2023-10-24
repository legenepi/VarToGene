#!/usr/bin/env Rscript

#Rationale: Run coloc and coloc.susie

suppressMessages(library(tidyverse))
library(coloc)
suppressMessages(library(Rfast))

args     = commandArgs(trailingOnly = TRUE)
tissue   = args[1]
cred_set = args[2]
gene     = args[3]

tmp_path="/scratch/gen1/nnp5/Var_to_Gen_tmp/"

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

############################
## Read asthma GWAS sumstat, eQTL GWAS, and LD matrix
############################
GWAS = read_delim(paste0(tmp_path, cred_set, "_GWASpairs.txt.gz")) %>%
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

eqtlGWAS = read_delim(paste0(tmp_path, cred_set, "_", tissue, "_pairs.txt.gz")) %>%
  mutate(ID      = gsub(x=ID, pattern=":", replacement="_"),
         varbeta = se^2) %>%
  arrange(chrom, pos)  %>%
  drop_na(beta) %>%
  filter(gene_id==gene, 
         ID %in% GWAS$snp) %>%
  rename(snp=ID, position=pos, MAF=maf)

GWAS %>% filter(snp %in% eqtlGWAS$snp) -> GWAS


############################
#Do colocalisation ONLY IF GWAS and eqtlGWAS (eQTL tissue-gene) contains pvalue <= 5x10-6.
############################
GWAS_sign <- GWAS %>% filter(pval <= 0.000005)
eqtlGWAS_sign <- eqtlGWAS %>% filter(pval <= 0.000005)
if (dim(eqtlGWAS_sign)[1] < 1){
    stop(paste0(signal," ",tissue, " ", gene, " GTExv8: No colocalisation is possible because eQTL data has pvalue >= 5x10-6."))
}
if (dim(GWAS_sign)[1] < 1){
    stop(paste0(signal," ",tissue, " ", gene, " GTExv8: No colocalisation is possible because GWAS data has pvalue >= 5x10-6."))
} else {
    paste0(signal," ",tissue, " ", gene, " GTExv8: Starting colocalisation because both GWAS and eQTL data has pvalue <= 5x10-6.")
}

############################
# * if also want coloc results
############################
coloc_all = coloc.abf(dataset1=list(beta=GWAS$beta, varbeta=GWAS$varbeta,
                                   N=GWAS$N, type="cc", MAF=GWAS$MAF, snp=GWAS$snp),
                     dataset2=list(beta=eqtlGWAS$beta, varbeta=eqtlGWAS$varbeta, 
                                   N=eqtlGWAS$N, type="quant", MAF=eqtlGWAS$MAF, snp=eqtlGWAS$snp))

saveRDS(coloc_all, paste0(tmp_path,"results/gtex/",cred_set, "_", tissue, "_", gene, "_all_coloc.rds"))

#look at the results:
coloc <- readRDS(paste0(tmp_path,"results/gtex/", cred_set, "_", tissue, "_", gene, "_all_coloc.rds"))
#sensitivity plot:
sensitivity(coloc,"H4 > 0.8")
#


############################
##COLOC.SUSIE to allow presence of multiple causal variants
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
    stop(paste0(signal," ",tissue, " ", gene, " GTExv8: No coloc.susie is possible because GWAS data has pvalue >= 5x10-6."))
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
  mutate(beta=ifelse(allele1!=ref, -(beta), beta)) # Check allele alignment

#Check eqtlGWAS has at least one variant with pval <= 0.000005:
eqtlGWAS_sign <- eqtlGWAS_df %>% filter(pval <= 0.000005)
if (dim(eqtlGWAS_sign)[1] < 1){
    stop(paste0(signal," ",tissue, " ", gene, " GTExv8: No coloc.susie is possible because eQTL-GWAS data has pvalue >= 5x10-6."))
} else {
    paste0(signal," ",tissue, " ", gene, " GTExv8: Starting coloc.susie because both GWAS and eQTL data has pvalue <= 5x10-6.")
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
  cat(paste0(cred_set, "_", tissue, "_", gene, ": ok", "\n"), 
      file=paste0(tmp_path,"GWAS_coloc_susie_check.txt"), append=TRUE)
} else {
  cat(paste0(cred_set, "_", tissue, "_", gene, ": warning", "\n"), 
      file=paste0(tmp_path,"GWAS_coloc_susie_check.txt"), append=TRUE)
}


############################
## Run SuSie
############################
susie_GWAS = runsusie(GWAS, r2.prune=0.2, check_R=FALSE)
susie_eQTLGWAS = runsusie(eqtlGWAS, r2.prune=0.2, check_R=FALSE)
susie_all = coloc.susie(susie_GWAS, susie_eQTLGWAS)

saveRDS(susie_all, paste0(tmp_path,"results/gtex/", cred_set, "_", tissue, "_", gene, "_all_susie.rds"))
write_tsv(susie_all$summary %>% as_tibble,
          paste0(tmp_path, "results/gtex/", cred_set, "_", tissue, "_", gene, "_susie_summary.tsv"))