#!/usr/bin/env Rscript

#Rationale: Variant Annotation with FAVOR for additional credset SNPs March 2025.

sink(stderr())

suppressMessages(library(tidyverse))
suppressMessages(library(data.table))


favor_file <- "/alice-home/3/n/nnp5/PhD/PhD_project/Var_to_Gene/input/Additional_credset_snps_March2025/20250310_FAVOR_output_additional_credset_SNPs_processed.csv.gz"
favor <- fread(favor_file,na.strings = c("",NA))

favor.digest <- favor %>% select("VariantVcf", "Chromosome", "Position", "Rsid", "GenecodeComprehensiveCategory",
                                 "GenecodeComprehensiveInfo", "GenecodeComprehensiveExonicCategory", "GenecodeComprehensiveExonicInfo",
                                 "CagePromoter", "CageEnhancer", "Genehancer", "SuperEnhancer", "CaddPhred", "FathmmXf",
                                 "ApcProteinFunctionV3", "ApcConservationV2", "ApcEpigeneticsActive", "ApcEpigeneticsRepressed",
                                 "ApcEpigeneticsTranscription", "ApcLocalNucleotideDiversityV3", "ApcMutationDensity",
                                 "ApcTranscriptionFactor", "ApcMappability", "Clnsig",
                                 "Clnsigincl", "Clndn", "Clndnincl", "Clnrevstat", "Origin",
                                 "Clndisdb", "Clndisdbincl", "Geneinfo")

colnames(favor.digest) <- make.names(colnames(favor.digest), unique=TRUE)

#FAVOR - Nearest Gene to variants with highest PIP in the locus:
#Identify whether variants cause protein coding changes using Gencode genes definition systems
#it will label the gene name of the variants has impact, if it is intergenic region
#the nearby gene name will be labeled in the annotation.


###Nearest genes:
#in Excel, I created the file with PIP-Sentinel for each credible per locus - if multiple credset, multiple PIP-Sentinel.
sentinel <- fread("input/Additional_credset_snps_March2025/Locus_PIPSentinel.txt")
sentinel <- sentinel %>% rename(Rsid=PIP_Sentinel)
favor.digest.NG <- favor.digest %>% select(Chromosome, Position, Rsid, "GenecodeComprehensiveInfo") %>% rename(Nearest_gene="GenecodeComprehensiveInfo")
sentinel_NG <- sentinel %>% left_join(favor.digest.NG, by = "Rsid") %>% rename(posb38=Position)
#for rs148639908: as the Rsid is not in favor, I manually added the gene:
favor.digest %>% filter(VariantVcf == "6-90253895-A-AT")
sentinel_NG$Chromosome[7] = 6
sentinel_NG$posb38[7] = as.numeric(90253895)
sentinel_NG$Nearest_gene[7] = "BACH2"
fwrite(sentinel_NG,"output/Additional_credset_snps_March2025_output/PIP_sentinels_nearestgenes_additionalcredset_March2025.txt",sep="\t",quote=F)
nearest_gene <- unique(unlist(strsplit(sentinel_NG$Nearest_gene,",")))

#Categorical annotation:
##Exploratory:
summary(as.factor(favor.digest$"GenecodeComprehensiveCategory"))
summary(as.factor(favor.digest$"GenecodeComprehensiveInfo"))
##Info for coding variants:
summary(as.factor(favor.digest$"GenecodeComprehensiveExonicCategory"))
summary(as.factor(favor.digest$"GenecodeComprehensiveExonicInfo"))

#FAVOR functional:
#FANTOM5 CAGE promoter or enhancer, Functional Integrative score > 15, ClinVar
#CAGE promoter or enhancer (from FANTOM5):
fantom5 <- favor.digest %>% filter(!is.na(CagePromoter) | !is.na(CageEnhancer))
fantom5_genes <- unique(fantom5$GenecodeComprehensiveInfo)
fantom5_df <- fantom5 %>% select(Chromosome, Position, GenecodeComprehensiveInfo, CagePromoter, CageEnhancer) %>% rename(Nearest_gene="GenecodeComprehensiveInfo")
fwrite(fantom5_df,"output/Additional_credset_snps_March2025_output/fnc_annot_fantom5_additionalcredset_March2025.txt",sep="\t",quote=F,row.names=F)


#Functional integrative score:
##Integrative score > 15:
inscores <- favor.digest %>% filter(CaddPhred > 15 | ApcProteinFunctionV3 > 15 | ApcConservationV2 > 15 | ApcEpigeneticsActive > 15 | ApcEpigeneticsRepressed > 15 |
ApcEpigeneticsTranscription > 15 | ApcLocalNucleotideDiversityV3 > 15 | ApcMutationDensity > 15 | ApcTranscriptionFactor > 15 | ApcMappability > 15)
inscores_genes <- unique(inscores$GenecodeComprehensiveInfo)

#Among variants with integrative score > 15 variants:
inscores.aPCs <- inscores %>% select("CaddPhred","ApcProteinFunctionV3", "ApcConservationV2", "ApcEpigeneticsActive", "ApcEpigeneticsRepressed",
                                 "ApcEpigeneticsTranscription", "ApcLocalNucleotideDiversityV3", "ApcMutationDensity",
                                 "ApcTranscriptionFactor", "ApcMappability")

inscores.aPCs_df <- inscores %>% select(Chromosome, Position, GenecodeComprehensiveInfo, "CaddPhred","ApcProteinFunctionV3", "ApcConservationV2", "ApcEpigeneticsActive", "ApcEpigeneticsRepressed",
                                 "ApcEpigeneticsTranscription", "ApcLocalNucleotideDiversityV3", "ApcMutationDensity",
                                 "ApcTranscriptionFactor", "ApcMappability")
fwrite(inscores.aPCs_df,"output/Additional_credset_snps_March2025_output/fnc_annot_inscores_additionalcredset_March2025.txt",sep="\t",quote=F,row.names=F)

#Clinical annotation:
clinvar <- favor.digest %>% filter(!is.na(Clnsig))
clin_genes <- unique(clinvar$Geneinfo)
clinvar_df <- clinvar %>% select(Chromosome, Position, Clnsig, Geneinfo)
fwrite(clinvar_df,"output/Additional_credset_snps_March2025_output/fnc_annot_clinvar_additionalcredset_March2025.txt",sep="\t",quote=F,row.names=F)

#Genes from Variant Annotation:
#fantom5_genes, inscores_genes, clin_genes
##Polish the genes (remove GeneID, remove ensembl, remove where mutation occurs, keep only gene name)
#fantom_genes: everything inside parenthesis needs to be deleted, and divide genes listed together
fantom5_genes <- gsub("\\(.*?\\)","",fantom5_genes) %>% unique()
fantom5_genes <- unlist(strsplit(fantom5_genes, "[\\s*,\\s*\\;]")) %>% unique()

#inscores_genes: everything inside parenthesis needs to be deleted, and divide genes listed together
inscores_genes <- gsub("\\(.*?\\)","",inscores_genes) %>% unique()
inscores_genes <- unlist(strsplit(inscores_genes, "\\s*,\\s*"))
inscores_genes <- unlist(strsplit(inscores_genes, ';', fixed=TRUE)) %>% unique()

#clin_genes: separate element if '|' present, and then remove GeneID-after the ':'
clin_genes <- unlist(strsplit(clin_genes, '|', fixed=TRUE))
clin_genes <- gsub(":.*", "", clin_genes) %>% unique()

varannot_genes <- as.data.frame(c(fantom5_genes,inscores_genes,clin_genes)) %>% unique()

#Save the genes:
fwrite(as.data.frame(nearest_gene),"input/Additional_credset_snps_March2025/nearest_genes_raw_additionalcredset_March2025.txt",row.names=FALSE, col.names=FALSE,quote=F)
fwrite(varannot_genes,"input/Additional_credset_snps_March2025/varannot_genes_additionalcredset_March2025.txt",row.names=FALSE, col.names=FALSE,quote=F)
#create /alice-home/3/n/nnp5/PhD/PhD_project/Var_to_Gene/input/Additional_credset_snps_March2025/var2genes_raw_additionalcredset_March2025.xlsx and add nearest_genes and varannot_genes sheets.