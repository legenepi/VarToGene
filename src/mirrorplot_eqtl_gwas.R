

#Ratioanle: in case I want to do mirrorplots, but I need more info from origianl summary stats.

library(mirrorplot)
df <- fread("SA_12_56435504_G_C_Skin_Not_Sun_Exposed_Suprapubic_lookup.txt.gz")
df <- df %>% rename(chr = chrom.gwas, pos = pos.gwas, trait1_p = pval.gwas, trait2_p = pval.eqtl, rsid = snpid)
mirrorplot(df, CHR = 12, START = 55935504, END = 56935504)
