#!/usr/bin/env Rscript

#Rationale: create table with results from each individual analysis:
#evidence	signal	gene	locus	start	end	chr	pos	trait	effect	other	eaf	Z	P	Novel	width

library(tidyverse)
library(data.table)
library(writexl)

#Cred_set variants:
#in credset, a1 is actually allele2 in gwas; a2 is actually allele1 in gwas
credset <- fread("/data/gen1/UKBiobank_500K/severe_asthma/Noemi_PhD/data/replsugg_valid_credset.txt") %>% select(-PIP_finemap,-PIP_susie) %>% rename(posb37=position,a1=allele2,a2=allele1,chr=chromosome)
#add posb38:
b38 <- fread("input/cs_vars_liftover_output.bed")
colnames(b38) <- c("chr_b38","posb38","pos1_b38")
credset_b38 <- cbind(credset,b38) %>% select(-chr_b38, -pos1_b38) %>%  relocate(posb38, .after = posb37)

#add gwas info to credset with also posb38:
gwas <- fread("/data/gen1/UKBiobank_500K/severe_asthma/Noemi_PhD/data/maf001_broad_pheno_1_5_ratio_betase_input_mungestat") %>% rename(chr=b37chr,posb37=bp)
credset_gwas <- credset_b38 %>% left_join(gwas,by=c("chr","posb37","a1","a2","snpid"))
credset_gwas <- credset_gwas %>% arrange(chr,posb37)
credset_gwas$sentinel <- paste0(credset_gwas$chr,"_",credset_gwas$posb37,"_",credset_gwas$a2,"_",credset_gwas$a1)


##Nearest gene:
ng <- fread("output/PIP_sentinels_nearestgenes") %>% select(locus,sentinel,Nearest_gene)
ng$evidence <- as.factor("nearest_gene")
credset_gwas_ng <- credset_gwas %>% left_join(ng, by=c("sentinel","locus"))
##extract only variants with evidence for nearest gene - the sentinel ones cus I did nearest gene by locus:
credset_gwas_ng2 <- credset_gwas_ng %>% filter(!is.na(evidence))

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

##Merge together the three functional criteria:
col_for_join <- c("locus","snpid","chr","posb37","posb38","a2","a1","PIP_average","LOG_ODDS","se","eaf","pval","MAF","sentinel","evidence")
fantom5_inscores <- full_join(credset_gwas_fantom5,credset_gwas_inscores.aPCs,by=col_for_join)
fantom5_inscores_clinvar <- full_join(fantom5_inscores,credset_gwas_clinvar,by=col_for_join) %>%
                            filter(!is.na(evidence)) %>%
                            unite(col='Gene', c('Nearest_gene', 'Genecode.Comprehensive.Info', 'Gene.Reported'),na.rm=TRUE)

##eQTL colocalisation:
#GTExV8:
gtex_coloc <- fread("output/coloc_asthma_GTEx.tsv") %>% select(snp, gene, tissue)
gtex_colocsusie <- fread("output/colocsusie_asthma_GTEx.tsv") %>% select(snp, gene, tissue)
gtex <- rbind(gtex_coloc, gtex_colocsusie)  %>% rename(sentinel_gtex=snp)
gtex$evidence <- as.factor("eqtl_gtex")
credset_gwas$sentinel_gtex <- paste0(credset_gwas$chr,"_",credset_gwas$posb37,"_",credset_gwas$a1,"_",credset_gwas$a2)
credset_gwas_gtex <- credset_gwas %>% left_join(gtex, by="sentinel_gtex")
credset_gwas_gtex2 <- credset_gwas_gtex %>% filter(!is.na(evidence))

#eqtlGen:
eqtlgen<- fread("output/coloc_asthma_eqtlgen.tsv") %>% select(snp, gene, tissue) %>% rename(sentinel_eqtlgen=snp)
eqtlgen$evidence <- as.factor("eqtl_eqtlgen")
credset_gwas$sentinel_eqtlgen <- paste0(credset_gwas$chr,"_",credset_gwas$posb37,"_",credset_gwas$a1,"_",credset_gwas$a2)
credset_gwas_eqtlgen <- credset_gwas %>% left_join(eqtlgen, by="sentinel_eqtlgen")
credset_gwas_eqtlgen2 <- credset_gwas_eqtlgen %>% filter(!is.na(evidence))

#UBCLung:
ubclung <- fread("output/coloc_asthma_ubclung.tsv") %>% select(snp, gene_id, tissue) %>% rename(sentinel_ubclung=snp)
ubclung$evidence <- as.factor("eqtl_ubclung")
credset_gwas$sentinel_ubclung <- paste0(credset_gwas$chr,"_",credset_gwas$posb37,"_",credset_gwas$a1,"_",credset_gwas$a2)
credset_gwas_ubclung <- credset_gwas %>% left_join(ubclung, by="sentinel_ubclung")
credset_gwas_ubclung2 <- credset_gwas_ubclung %>% filter(!is.na(evidence))

##pQTL look-up:
#UKB pQTL:
ukbpqtl <- fread("output/lookup_ukbpqtl.txt")
ukbpqtl <- ukbpqtl %>% separate(LOCUS, c("locus", "sentinel_ukbpqtl"), sep = "/")
ukbpqtl <- ukbpqtl %>% rename(chr=CHROM,posb38=GENPOS)
ukbpqtl$evidence <- as.factor("pqtl_ukb")
credset_gwas_ukbpqtl <- credset_gwas %>% left_join(ukbpqtl, by=c("locus","chr","posb38"))
credset_gwas_ukbpqtl2 <- credset_gwas_ukbpqtl %>% filter(!is.na(evidence))

#SCALLOP:
#Colnmaes:#Chrom	Pos	MarkerName	Allele1	Allele2	Freq1	FreqSE	Effect	StdErr	P-value	Direction	TotalSampleSize Gene
scallop <- fread("output/scallop_ukbpqtl.txt",head=F,sep="\t")
scallop <- scallop %>% separate(V12, c("TotalSampleSize","Gene"),sep=" ")
colnames(scallop) <- c("Chrom","Pos","MarkerName","Allele1","Allele2","Freq1","FreqSE","Effect", "StdErr","P-value","Direction","TotalSampleSize","Gene")
scallop <- scallop %>% rename(chr=Chrom,posb37=Pos)
scallop$evidence <- as.factor("pqtl_scallop")
credset_gwas_scallop <- credset_gwas %>% left_join(scallop, by=c("chr","posb37"))
credset_gwas_scallop2 <- credset_gwas_scallop %>% filter(!is.na(evidence))

##Mouse_ko
mko <- fread("input/mko_results_500kb.csv") %>% filter(overlap == TRUE) %>% select(-Position, -overlap, -MKO) %>% rename(posb37=pos) %>% arrange(chr,posb37)
mko$evidence <- as.factor("mouse_ko")
credset_gwas_mko <- credset_gwas %>% left_join(mko, by=c("locus","sentinel","chr","posb37"))
credset_gwas_mko2 <- credset_gwas_mko %>% filter(!is.na(evidence))

##Rare disease
raredis <- fread("input/raredis_results_500kb_by_gene.csv")
colnames(raredis) <- c("Symbol","sentinel","frequency_disease","distance","name_disease","HPOterm")
raredis$evidence <- as.factor("rare_disease")
credset_gwas_raredis <- credset_gwas %>% left_join(raredis, by= "sentinel")
credset_gwas_raredis2 <- credset_gwas_raredis %>% filter(!is.na(evidence))

#Save all results for each analysis into a unique xlsx file:
df_list <- list(credset_gwas_ng2,fantom5_inscores_clinvar,credset_gwas_gtex2,credset_gwas_eqtlgen2,
credset_gwas_ubclung2,credset_gwas_ukbpqtl2,credset_gwas_scallop2,credset_gwas_raredis2,credset_gwas_mko2)
write_xlsx(df_list,path = "src/report/var2gene_full.xlsx", col_names = TRUE, format_headers = TRUE)

##Rare variant ExWAS
##Single rare-variant:

##Gene-collpasing rare variant:

