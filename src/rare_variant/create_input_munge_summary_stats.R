#!/usr/bin/env Rscriptcd all

#Rationale filter out rare variants form ExWAS.

library(data.table)
library(tidyverse)

args = commandArgs(trailingOnly=TRUE)
df_file <- args[1]
pheno <- as.character(args[2])

df <- fread(df_file,header=F,sep="\t")

#"#chrom"          "pos"             "rsid"            "ref"
#"alt"             "neg_log_pvalue"  "beta"            "stderr_beta"
#"alt_allele_freq"
#CHROM GENPOS ID ALLELE0 ALLELE1 A1FREQ INFO N TEST BETA SE CHISQ LOG10P
#REGENIE: reference allele (allele 0), alternative allele (allele 1)
#REGENIE:  estimated effect sizes (for allele 1 on the original scale)

colnames(df) <- c("CHROM", "GENPOS", "ID", "ALLELE0", "ALLELE1", "LOG10P", "BETA", "SE", "A1FREQ")
dfclean <- df %>% select(CHROM,ID,GENPOS,ALLELE1,ALLELE0,A1FREQ,LOG10P,BETA,SE)

dfclean$pval <- 10^(-as.numeric(dfclean$LOG10P))

#select column in order:
dfclean <- dfclean %>% select(ID,CHROM,GENPOS,ALLELE1,ALLELE0,BETA,SE,A1FREQ,pval)

#rename columns :
colnames(dfclean) <- c("snpid", "chr", "bp38", "a1", "a2", "LOG_ODDS", "se", "eaf", "pval")

#make sure eaf column is numeric:
dfclean$eaf <- as.numeric(dfclean$eaf)

#create MAF columns and filtered out for MAF >= 0.01:
dfclean <- dfclean %>% mutate(MAF = ifelse(dfclean$eaf <= 0.5, dfclean$eaf, 1 - dfclean$eaf))
dfclean_rarevar <- dfclean %>% filter(MAF < 0.01)

#save gwas results for variants with MAF < 0.01:
write.table(dfclean_rarevar,paste0("input/rare_variant/single_rarevar_EwWAS/",pheno,"rarevar_betase_input_mungestat"),quote=FALSE,sep="\t",row.names=FALSE)

#save gwas results for variants with MAF < 0.01 and pvalue <= 0.000005:
dfclean_rarevar_sugg <- dfclean_rarevar %>% filter(pval <= 0.000005)
write.table(dfclean_rarevar_sugg,paste0("input/rare_variant/single_rarevar_EwWAS/",pheno,"rarevar_suggestive"),quote=FALSE,sep="\t",row.names=FALSE)