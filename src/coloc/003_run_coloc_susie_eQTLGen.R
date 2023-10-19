suppressMessages(library(tidyverse))
library(coloc)
suppressMessages(library(Rfast))

args     = commandArgs(trailingOnly = TRUE)
cred_set = args[1]
gene     = args[2]

cred_set %>% 
  as_tibble %>% 
  separate(value, c("credible", "anc"), sep="_") %>% 
  unite("credible", credible:anc) %>% 
  pull -> credible_set

cred_set %>% 
  as_tibble %>% 
  separate(value, c("credible", "anc", "chr", "pos", "a1", "a2"), sep="_") %>% 
  unite("signal", chr:a2) %>% 
  pull -> signal

############################
## Read neuropathic pain and eQTL GWAS, and LD matrix 
############################
npGWAS = read_delim(paste0("/data/gen1/np_gwas/np_GWAS_3/analysis/fine_mapping/", cred_set, ".txt")) %>% 
  select(-snp) %>%
  arrange(chr, pos) %>%
  drop_na(beta) %>%
  mutate(varbeta        = SE^2, 
         MAF            = if_else(AF.alt <0.5, AF.alt, 1-AF.alt), 
         new_variant_id = paste0(chr, "_", pos, "_", pmin(ref, alt), "_", pmax(ref, alt))) %>%
  rename(snp=new_variant_id, position=pos, N=num) %>%
  distinct(snp, .keep_all=TRUE)

eqtlGWAS = read_delim(paste0("/scratch/gen1/atw20/pain/results/coloc/eqtlgen/", cred_set, "_pairs.txt.gz")) %>% 
  mutate(ID      = gsub(x=ID, pattern=":", replacement="_"), 
         varbeta = se^2, 
         MAF     = if_else(AssessedAllele_freq <0.5, AssessedAllele_freq, 1-AssessedAllele_freq)) %>% 
  arrange(chrom, pos)  %>%
  drop_na(beta) %>%
  filter(GeneSymbol==gene, 
         ID %in% npGWAS$snp) %>%
  rename(snp=ID, position=pos, N=NrSamples)

npGWAS %>% filter(snp %in% eqtlGWAS$snp) -> npGWAS

############################
# * if also want coloc results
############################
coloc_all = coloc.abf(dataset1=list(beta=npGWAS$beta, varbeta=npGWAS$varbeta, 
                                    N=npGWAS$N, type="cc", MAF=npGWAS$MAF, snp=npGWAS$snp), 
                      dataset2=list(beta=eqtlGWAS$beta, varbeta=eqtlGWAS$varbeta, 
                                    N=eqtlGWAS$N, type="quant", MAF=eqtlGWAS$MAF, snp=eqtlGWAS$snp))

saveRDS(coloc_all, paste0("/scratch/gen1/atw20/pain/results/coloc/eqtlgen/coloc_results/", cred_set, "_", gene, "_all_coloc.rds"))

############################
# 1) Format LD matrix
############################
ld = data.table::fread(
  paste0("/scratch/gen1/atw20/pain/results/coloc/eqtlgen/", signal, "_", credible_set, "_whole_blood.raw")
) %>% 
  select(!c(FID:PHENOTYPE))

ld %>% 
  names %>% 
  as_tibble %>% 
  separate(value,c("c", "p", "a1", "a2")) %>% 
  unite("snp", c:a2, sep="_") %>% 
  pull -> snps

ld = t(ld) %>%
  as_tibble %>%
  mutate(var=apply(., 1, var, na.rm=TRUE),
         snp=snps) %>%
  filter(var!=0) %>%
  filter(snp %in% npGWAS$snp) %>%
  filter(snp %in% eqtlGWAS$snp)

N = ncol(ld)-2 
X = t(as.matrix(ld[, 1:N]))
colnames(X) = ld %>% 
  select(snp) %>% 
  pull # Sdd names to facilitate filtering NAs (still some NAs even with pairwise.complete.obs)
LDmatrix = cor(X, use="pairwise.complete.obs") # Pearson's correlation, takes some time...
LDmatrix[is.na(LDmatrix)] <- 0 

############################
# 2) Format neuropathic GWAS
############################
bim = read_tsv(
  paste0("/scratch/gen1/atw20/pain/results/coloc/eqtlgen/", signal, "_", credible_set, "_whole_blood.bim"), 
  col_names=c("chr", "snp", "morgans", "position", "allele1", "allele2")
) %>% 
  select(!morgans) # Read bim file to check same set of SNPs and allele alignment

npGWAS %>%  
  filter(snp %in% colnames(LDmatrix)) %>% # Must be the same set of SNPs
  inner_join(bim) %>%
  mutate(beta=ifelse(allele1!=ref, -(beta), beta)) %>% # Check allele alignment
  as.list -> npGWAS

npGWAS$type = "cc"
npGWAS$LD = LDmatrix
N = as.integer(mean(npGWAS$N))
npGWAS$N = N

############################
# 3) Format eQTL GWAS  
############################
eqtlGWAS  %>% 
  filter(snp %in% colnames(LDmatrix)) %>% # Must be the same set of SNPs 
  inner_join(bim) %>%
  mutate(beta=ifelse(allele1!=AssessedAllele, -(beta), beta)) %>% # Check allele alignment
  as.list -> eqtlGWAS

eqtlGWAS$type = "quant"
eqtlGWAS$LD = LDmatrix
N = as.integer(mean(eqtlGWAS$N))
eqtlGWAS$N = N

############################
# Check datasets are ok
############################
check_dataset(eqtlGWAS, req="LD") -> eqtlGWAS_check
check_dataset(npGWAS, req="LD") -> npGWAS_check

if (is.null(npGWAS_check) & is.null(eqtlGWAS_check)){
  cat(paste0(cred_set, "_", tissue, "_", gene, ": ok", "\n"), 
      file="/scratch/gen1/atw20/pain/results/coloc/eqtlgen/coloc_results/GWAS_coloc_susie_check.txt", append=TRUE)
} else {
  cat(paste0(cred_set, "_", tissue, "_", gene, ": warning", "\n"), 
      file="/scratch/gen1/atw20/pain/results/coloc/eqtlgen/coloc_results/GWAS_coloc_susie_check.txt", append=TRUE)
}

############################
## Run SuSie
############################
susie_npGWAS = runsusie(npGWAS, r2.prune=0.2, check_R=FALSE)
susie_eQTLGWAS = runsusie(eqtlGWAS, r2.prune=0.2, check_R=FALSE)
susie_all = coloc.susie(susie_npGWAS, susie_eQTLGWAS)

write_tsv(susie_all$summary, paste0("/scratch/gen1/atw20/pain/results/coloc/eqtlgen/coloc_results/", cred_set, "_", gene, "_susie_summary.tsv"))
saveRDS(susie_all, paste0("/scratch/gen1/atw20/pain/results/coloc/eqtlgen/coloc_results/", cred_set, "_", gene, "_all_susie.rds"))
write_tsv(susie_all$summary %>% as_tibble, 
          paste0("/scratch/gen1/atw20/pain/results/coloc/eqtlgen/coloc_results/", cred_set, "_", gene, "_susie_summary.tsv"))
