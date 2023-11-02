#!/usr/bin/env Rscript

#Rationale: Create GTExv8 eQTL files to use for colocalisation, filter genes for each eQTL files, find GWAS sumstat lookup
#in GTEXv8 for regions that I cannot do colocalisation, e.i. HLA regions.

library(data.table)
library(tidyverse)

args     = commandArgs(trailingOnly = TRUE)
tissue   = args[1]
cred_set = args[2]

tmp      = unlist(strsplit(cred_set, split="_"))
chr      = as.numeric(tmp[2])
location = as.numeric(tmp[3])

tmp_path="/scratch/gen1/nnp5/Var_to_Gen_tmp/"

##allele are place in alphabetical order ! (no dependency on effect allele, reference/alternative)

#eqtl:
#'Colon_Transverse' 'Colon_Sigmoid' 'Skin_Sun_Exposed_Lower_leg' 'Skin_Not_Sun_Exposed_Suprapubic' files are in my scratch;
#make an exeption for these files:
my_tissues <- c("Colon_Transverse", "Colon_Sigmoid", "Skin_Sun_Exposed_Lower_leg", "Skin_Not_Sun_Exposed_Suprapubic")
if (tissue %in% my_tissues ){
    eqtl = fread(paste0(tmp_path,"liftover_gtexv8/", tissue, ".v8.EUR.allpairs.chr", chr, ".hg19.txt.gz"),
             select=c("phenotype_id", "CHROM", "POS", "REF", "ALT", "maf", "ma_samples", "slope", "slope_se", "pval_nominal"),
             col.names=c("gene_id", "chrom", "pos", "ref", "alt", "maf", "N", "beta", "se", "pval"),
             showProgress=FALSE) %>%
        mutate(ID = paste0(chrom, ":", pos, "_", pmin(ref, alt), "_", pmax(ref, alt))) %>%
        relocate(ID, .before="gene_id") %>%
        filter(between(pos, location-5e5, location+5e5), chrom==chr)
    } else {
eqtl = fread(paste0("/data/gen1/ACEI/colocalisation_datasets/eQTL/GTeX/", tissue, ".v8.EUR.allpairs.chr", chr, ".hg19.txt.gz"), 
             select=c("phenotype_id", "CHROM", "POS", "REF", "ALT", "maf", "ma_samples", "slope", "slope_se", "pval_nominal"), 
             col.names=c("gene_id", "chrom", "pos", "ref", "alt", "maf", "N", "beta", "se", "pval"), 
             showProgress=FALSE) %>%
  mutate(ID = paste0(chrom, ":", pos, "_", pmin(ref, alt), "_", pmax(ref, alt))) %>%
  relocate(ID, .before="gene_id") %>%
  filter(between(pos, location-5e5, location+5e5), chrom==chr)
}

genes = eqtl %>% distinct(gene_id)

write_delim(x=eqtl, file=paste0(tmp_path, cred_set, "_", tissue, "_pairs.txt.gz"))
write_delim(x=genes, file=paste0(tmp_path, cred_set, "_", tissue, "_genes.txt"), col_names=FALSE)

#import gwas sumstat:
gwas_sumstat <- fread("/data/gen1/UKBiobank_500K/severe_asthma/Noemi_PhD/data/maf001_broad_pheno_1_5_ratio_betase_input_mungestat")
gwas_sumstat <- gwas_sumstat %>% setnames(c("b37chr","bp","a1","a2"),c("chrom","position","allele1","allele2"))

#credible set:
cs <- fread(paste0("/scratch/gen1/nnp5/Var_to_Gen_tmp/", cred_set))
cs <- cs %>% setnames("chromosome","chrom")

#merge credible set with gwas sumstat:
cs_gwas_sumstat <- inner_join(cs,gwas_sumstat, by=c("chrom","position","allele1","allele2","snpid"))

#merge credible set-gwas sumstat with eqtl:
cs_gwas_sumstat_eqtl <- cs_gwas_sumstat %>% setnames(c("position","LOG_ODDS"),c("pos","beta")) %>%
  mutate(ID = paste0(chrom, ":", pos, "_", pmin(allele1, allele2), "_", pmax(allele1, allele2))) %>%
  relocate(ID, .before="chrom") %>%
  distinct() %>%
  left_join(eqtl, by="ID", suffix=c(".gwas", ".eqtl")) %>%
  mutate(beta.eqtl = if_else(alt == allele1, beta.eqtl, -beta.eqtl))

write_delim(x=cs_gwas_sumstat_eqtl, file=paste0(tmp_path, cred_set, "_", tissue, "_lookup.txt.gz"))