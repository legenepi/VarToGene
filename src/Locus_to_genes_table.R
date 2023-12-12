#!/usr/bin/env Rscript

#Rationale: create table with results from each individual analysis:
#evidence	signal	gene	locus	start	end	chr	pos	trait	effect	other	eaf	Z	P	Novel	width

library(tidyverse)
library(data.table)

#Cred_set variants:
#in credset, a1 is actually allele2 in gwas; a2 is actually allele1 in gwas
credset <- fread("/data/gen1/UKBiobank_500K/severe_asthma/Noemi_PhD/data/replsugg_valid_credset.txt") %>% select(-PIP_finemap,-PIP_susie) %>% rename(posb37=position,a1=allele2,a2=allele1,chr=chromosome)
#creset posb38:
b38 <- fread("input/cs_vars_liftover_output.bed")
colnames(b38) <- c("chr_b38","posb38","pos1_b38")
credset_b38 <- cbind(credset,b38) %>% select(-chr_b38, -pos1_b38)

gwas <- fread("/data/gen1/UKBiobank_500K/severe_asthma/Noemi_PhD/data/maf001_broad_pheno_1_5_ratio_betase_input_mungestat") %>% rename(chr=b37chr,posb37=bp)
credset_gwas <- credset_b38 %>% left_join(gwas,by=c("chr","posb37","a1","a2","snpid"))
credset_gwas <- credset_gwas %>% arrange(chr,posb37)
credset_gwas$sentinel <- paste0(credset_gwas$chr,"_",credset_gwas$posb37,"_",credset_gwas$a2,"_",credset_gwas$a1)


##Nearest gene:
ng <- fread("output/PIP_sentinels_nearestgenes") %>% select(locus,sentinel,Nearest_gene)
ng$evidence <- as.factor("nearest_gene")
credset_gwas_ng <- credset_gwas %>% left_join(ng, by=c("sentinel","locus"))

##Functional annotation:
#fantom5
fantom5 <- fread("output/fnc_annot_fantom5") %>% rename(chr=Chromosome,posb38=Position)
fantom5$evidence <- as.factor("fantom5")
credset_gwas_fantom5 <- credset_gwas %>% left_join(fantom5, by=c("chr","posb38"))

#integrative scores:
inscores.aPCs <- fread("output/fnc_annot_inscores") %>% rename(chr=Chromosome,posb38=Position)
inscores.aPCs$evidence <- as.factor("inscores.aPCs")
credset_gwas_inscores.aPCs <- credset_gwas %>% left_join(inscores.aPCs, by=c("chr","posb38"))

#ClinVar:
clinvar <- fread("output/fnc_annot_clinvar") %>% rename(chr=Chromosome,posb38=Position)
clinvar$evidence <- as.factor("clinvar")
credset_gwas_clinvar <- credset_gwas %>% left_join(clinvar, by=c("chr","posb38"))


##eQTL colocalisation:
#GTExV8:
gtex_coloc <- fread("output/coloc_asthma_GTEx.tsv") %>% select(snp, gene, tissue)
gtex_colocsusie <- fread("output/colocsusie_asthma_GTEx.tsv") %>% select(snp, gene, tissue)
gtex <- rbind(gtex_coloc, gtex_colocsusie)  %>% rename(sentinel_gtex=snp)
gtex$evidence <- as.factor("eqtl_gtex")
credset_gwas$sentinel_gtex <- paste0(credset_gwas$chr,"_",credset_gwas$posb37,"_",credset_gwas$a1,"_",credset_gwas$a2)
credset_gwas_gtex <- credset_gwas %>% left_join(gtex, by="sentinel_gtex")

#eqtlGen:
eqtlgen<- fread("output/coloc_asthma_eqtlgen.tsv") %>% select(snp, gene, tissue) %>% rename(sentinel_eqtlgen=snp)
eqtlgen$evidence <- as.factor("eqtl_eqtlgen")
credset_gwas$sentinel_eqtlgen <- paste0(credset_gwas$chr,"_",credset_gwas$posb37,"_",credset_gwas$a1,"_",credset_gwas$a2)
credset_gwas_eqtlgen <- credset_gwas %>% left_join(eqtlgen, by="sentinel_eqtlgen")

#UBCLung:
ubclung <- fread("output/coloc_asthma_ubclung.tsv") %>% select(snp, gene_id, tissue) %>% rename(sentinel_ubclung=snp)
ubclung$evidence <- as.factor("eqtl_ubclung")
credset_gwas$sentinel_ubclung <- paste0(credset_gwas$chr,"_",credset_gwas$posb37,"_",credset_gwas$a1,"_",credset_gwas$a2)
credset_gwas_ubclung <- credset_gwas %>% left_join(ubclung, by="sentinel_ubclung")

##pQTL look-up:
#UKB pQTL:
ukbpqtl <- fread("output/lookup_ukbpqtl.txt")
ukbpqtl <- ukbpqtl %>% separate(LOCUS, c("locus", "sentinel_ukbpqtl"), sep = "/")
ukbpqtl <- ukbpqtl %>% rename(chr=CHROM,posb38=GENPOS)
ukbpqtl$evidence <- as.factor("pqtl_ukb")
credset_gwas_ukbpqtl <- credset_gwas %>% left_join(ukbpqtl, by=c("locus","chr","posb38"))

#SCALLOP:
#Colnmaes:#Chrom	Pos	MarkerName	Allele1	Allele2	Freq1	FreqSE	Effect	StdErr	P-value	Direction	TotalSampleSize Gene
scallop <- fread("output/scallop_ukbpqtl.txt",head=F,sep="\t")
scallop <- scallop %>% separate(V12, c("TotalSampleSize","Gene"),sep=" ")
colnames(scallop) <- c("Chrom","Pos","MarkerName","Allele1","Allele2","Freq1","FreqSE","Effect", "StdErr","P-value","Direction","TotalSampleSize","Gene")
scallop <- scallop %>% rename(chr=Chrom,posb37=Pos)
credset_gwas_scallop <- credset_gwas %>% left_join(scallop, by=c("chr","posb37"))

##Mouse_ko


##Rare disease


##Rare variant ExWAS
##Single rare-variant:

##Gene-collpasing rare variant:

