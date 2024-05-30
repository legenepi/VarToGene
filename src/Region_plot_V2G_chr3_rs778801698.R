#!/usr/bin/env Rscript


#Rationale: plot of Genomic Locus and variant-to-gene mapping genes:
##variants in credible set, with lead variant as the one with highest PIP
##functional annotation for only credible set variants, using as by FAVOR
##R2 between variants in the region
##Gene highlighted by the variant-to-gene mapping analysis with the type of evidence

library(tidyverse)
library(data.table)
library(dplyr)

annot <- read.table("/data/gen1/UKBiobank_500K/severe_asthma/Noemi_PhD/data/FAVOR_credset_chrpos38_2024_05_14_chr3_49524027_50524027_rs778801698.txt.csv", sep=",", stringsAsFactors = FALSE, fill=TRUE, header=TRUE)
annot <- annot %>% select("VariantVcf","Chromosome","Position","GenecodeComprehensiveCategory") %>%
    rename(id="VariantVcf",chromosome="Chromosome",posb38="Position")
annot <- as.data.frame(sapply(annot, function(x) gsub("\"", "", x)))
annot <- annot %>% rename(Functional_annotation=GenecodeComprehensiveCategory)

#liftover online:
#awk -F "_" '{print "chr"$1,$2,$2}' /scratch/gen1/nnp5/Var_to_Gen_tmp/regionalplot/cs_variants_EUR_UKB_${region} > /home/n/nnp5/PhD/PhD_project/Var_to_Gene/input/cs_variants_EUR_UKB_lifover_input_${region}
credset <- fread("/scratch/gen1/nnp5/Var_to_Gen_tmp/regionalplot/cs_variants_EUR_UKB_chr3_49524027_50524027_rs778801698",header=FALSE)
credset <- credset %>% separate(V1, "_", into=c("chromosome", "posb37", "allele1", "allele2"))
b38 <- fread("/home/n/nnp5/PhD/PhD_project/Var_to_Gene/input/cs_variants_EUR_UKB_liftover_output_chr3_49524027_50524027_rs778801698.bed",header=FALSE)
credset_b38 <- cbind(credset,b38)
credset_b38 <- credset_b38 %>% rename(posb38=V2)
credset_b38$chromosome <- as.numeric(credset_b38$chromosome)
annot$chromosome <- as.numeric(annot$chromosome)
credset_b38$posb38 <- as.numeric(credset_b38$posb38)
annot$posb38 <- as.numeric(annot$posb38)
credset_annot <- inner_join(credset_b38, annot, by=c("chromosome","posb38"))
credset_annot$credset <- "1"
credset_annot$posb37 <- as.numeric(credset_annot$posb37)


gwas <- fread("/data/gen1/UKBiobank_500K/severe_asthma/Noemi_PhD/data/maf001_broad_pheno_1_5_ratio_betase_input_mungestat") %>% rename(chromosome=b37chr,posb37=bp)
gwas_credset_annot <- left_join(gwas,credset_annot,by=c("chromosome","posb37"))
gwas_credset_annot <- gwas_credset_annot %>% select(chromosome,posb37,pval,a1,a2,Functional_annotation,credset,snpid)
gwas_credset_annot <- gwas_credset_annot %>% mutate(Functional_annotation=ifelse(is.na(gwas_credset_annot$Functional_annotation), "unknown",gwas_credset_annot$Functional_annotation))
gwas_credset_annot$Functional_annotation <- as.factor(gwas_credset_annot$Functional_annotation)
gwas_credset_annot <- gwas_credset_annot %>% mutate(credset=ifelse(is.na(gwas_credset_annot$credset), 0, gwas_credset_annot$credset))
fwrite(gwas_credset_annot,"/scratch/gen1/nnp5/Var_to_Gen_tmp/regionalplot/gwas_credset_annot_chr3_49524027_50524027_rs778801698",quote=F,sep="\t",row.names=F,col.names=T)