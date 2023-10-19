library(data.table)
library(tidyverse)

args     = commandArgs(trailingOnly = TRUE)
cred_set = args[1]

tmp      = unlist(strsplit(cred_set, split="_"))
chr      = as.numeric(tmp[3])
location = as.numeric(tmp[4])
pheno    = paste0(tmp[1], "_", tmp[2])

eqtl = fread("/scratch/gen1/atw20/pain/results/coloc/eqtlgen/whole_blood_cis_eqtls_withAF.txt.gz", 
             showProgress=FALSE) %>%
  setnames(c("SNPChr", "SNPPos", "Zscore", "Pvalue"), c("chrom", "pos", "Z_score", "pval")) %>%
  mutate(ID = paste0(chrom, ":", pos, "_", pmin(AssessedAllele, OtherAllele), "_", pmax(AssessedAllele, OtherAllele))) %>%
  relocate(ID, .before="SNP") %>%
  filter(between(pos, location-5e5, location+5e5), chrom==chr)

genes = eqtl %>% distinct(GeneSymbol)

write_delim(x=eqtl, file=paste0("/scratch/gen1/atw20/pain/results/coloc/eqtlgen/", cred_set, "_pairs.txt.gz"))
write_delim(x=genes, file=paste0("/scratch/gen1/atw20/pain/results/coloc/eqtlgen/", cred_set, "_genes.txt"), col_names=FALSE)

cs = fread(paste0("/data/gen1/np_gwas/np_GWAS_3/analysis/fine_mapping/", cred_set, ".txt"), 
           drop=c("CHR_A", "BP_A", "CHR_B", "BP_B", "SNP_B", "id", "R2"), 
           showProgress=FALSE) %>% 
  setnames(c("chr", "SE", "P"), c("chrom", "se", "pval")) %>%
  mutate(ID = paste0(chrom, ":", pos, "_", pmin(ref, alt), "_", pmax(ref, alt))) %>%
  relocate(ID, .before="chrom") %>%
  distinct() %>%
  left_join(eqtl, by="ID", suffix=c(".gwas", ".eqtl")) %>%
  mutate(beta.eqtl    = if_else(alt == AssessedAllele, beta.eqtl, -beta.eqtl), 
         credible_set = pheno) %>%
  relocate(credible_set, .before="SNP_A")

write_delim(x=cs, file=paste0("/scratch/gen1/atw20/pain/results/coloc/eqtlgen/", cred_set, "_lookup.txt.gz"))


# P2200_AFR_21_20113911_C_T
# P2200.3_AFR_21_20113911_C_T
# P2200_EUR_11_112917378_C_T
# P2200.3_EUR_11_112917378_C_T
# P2200.1_EUR_6_32460831_G_T
# P2200.39_EUR_1_192923542_A_G
# P2200.40_EUR_5_113421943_A_G
# P2200.51_EUR_10_73933699_A_G
# P2200.511_EUR_1_64613556_A_G
