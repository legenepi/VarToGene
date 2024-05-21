#!/usr/bin/env Rscript

#Rationale: Variant Annotation with FAVOR - need a new format because the colnames have changed as FAVOR got updated

sink(stderr())

suppressMessages(library(tidyverse))
suppressMessages(library(data.table))
suppressMessages(library(corrplot))
suppressMessages(library(VennDiagram))
library(RColorBrewer)

args <- commandArgs(T)
favor_file <- args[1]
region <- args[2]
#region <- "chr3_49524027_50524027_rs77880169"
#favor_file <- "/alice-home/3/n/nnp5/PhD/PhD_project/Var_to_Gene/input/FAVOR_credset_chrpos38_2024_05_14_chr3_49524027_50524027_rs778801698.txt.csv"
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
##Retrieve Sentinel SNPs - highest PIP in the fine-mapped loci:
#cp /scratch/gen1/nnp5/Var_to_Gen_tmp/ukb_pqtl/cs_vars_liftover_output.bed /home/n/nnp5/PhD/PhD_project/Var_to_Gene/input/
sig_list_tmp <- fread("/data/gen1/UKBiobank_500K/severe_asthma/Noemi_PhD/data/replsugg_valid_credset_chr3_noMHC.txt")
setnames(sig_list_tmp,"chromosome","chr")
setnames(sig_list_tmp,"position","posb37")
sig_list_tmp$sentinel <- paste0(sig_list_tmp$chr,"_",sig_list_tmp$posb37,"_",sig_list_tmp$allele1,"_",sig_list_tmp$allele2)
locus <- unique(sig_list_tmp$locus)

sentinels <- data.frame(matrix(ncol = 4,nrow = 0))
colnames(sentinels) <- c("locus","sentinel","chr","posb37")
locus_sig_list <- sig_list_tmp %>% filter(locus == "3_rs778801698_49524027_50524027")
locus_sig_list <- locus_sig_list %>% filter(PIP_average == max(locus_sig_list$PIP_average)) %>% select(locus,sentinel,chr,posb37)
sentinels <- rbind(sentinels,locus_sig_list)

#liftover fo rhte senitnel to find b38:
posb38 <- fread("/home/n/nnp5/PhD/PhD_project/Var_to_Gene/input/hglft_genome_credset_vars_chr3_49524027_50524027_rs778801698.bed") %>% rename(posb38=V2) %>% select(posb38)
posb37 <- fread("/home/n/nnp5/PhD/PhD_project/Var_to_Gene/input/replsugg_valid_credset_input_liftover_b37_chr3") %>% rename(posb37=V2) %>% select(posb37)
pos <- cbind(posb38,posb37)
sentinel_b38 <- sentinels %>% left_join(pos, by= "posb37") %>% rename(Chromosome=chr, Position=posb38)


###Nearest genes:
favor.digest.NG <- favor.digest %>% select(Chromosome, Position, "GenecodeComprehensiveInfo") %>% rename(Nearest_gene="GenecodeComprehensiveInfo")
sentinel_b38_NG <- sentinel_b38 %>% left_join(favor.digest.NG, by = c("Chromosome", "Position")) %>% rename(posb38=Position)
fwrite(sentinel_b38_NG,paste0("output/PIP_sentinels_nearestgenes",region),sep="\t",quote=F)
nearest_gene <- unique(unlist(strsplit(sentinel_b38_NG$Nearest_gene,",")))

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
fwrite(fantom5_df,paste0("output/fnc_annot_fantom5",region),sep="\t",quote=F,row.names=F)


#Functional integrative score:
#Among all the variants, there are a correlations between aPC.Epigenetics.Active and aPC.Transcription.Factor (0.75), between aPC.Mutation.Density
#and aPC.Local.Nucleotide.Diversity (0.80); between CADD.phred and aPCs.Conservation (0.83).
favor.aPCs <- favor.digest %>% select("CaddPhred","ApcProteinFunctionV3", "ApcConservationV2", "ApcEpigeneticsActive", "ApcEpigeneticsRepressed",
                                 "ApcEpigeneticsTranscription", "ApcLocalNucleotideDiversityV3", "ApcMutationDensity",
                                 "ApcTranscriptionFactor", "ApcMappability")

favor.aPCs <- na.omit(favor.aPCs)
M <- cor(favor.aPCs)
head(round(M,3))
corrplot(M, method="number")

##Integrative score > 15:
inscores <- favor.digest %>% filter(CaddPhred > 15 | ApcProteinFunctionV3 > 15 | ApcConservationV2 > 15 | ApcEpigeneticsActive > 15 | ApcEpigeneticsRepressed > 15 |
ApcEpigeneticsTranscription > 15 | ApcLocalNucleotideDiversityV3 > 15 | ApcMutationDensity > 15 | ApcTranscriptionFactor > 15 | ApcMappability > 15)
inscores_genes <- unique(inscores$GenecodeComprehensiveInfo)

#Among variants with integrative score > 15 variants:
inscores.aPCs <- inscores %>% select("CaddPhred","ApcProteinFunctionV3", "ApcConservationV2", "ApcEpigeneticsActive", "ApcEpigeneticsRepressed",
                                 "ApcEpigeneticsTranscription", "ApcLocalNucleotideDiversityV3", "ApcMutationDensity",
                                 "ApcTranscriptionFactor", "ApcMappability")
M2 <- cor(inscores.aPCs)
head(round(M,3))
#Save the corplot x the report:
png(file = paste0("/home/n/nnp5/PhD/PhD_project/Var_to_Gene/output/corrplot_integrative_score_15", region, ".png"))
corplot_integrativescore <- corrplot(M2, method="number")
dev.off()

inscores.aPCs_df <- inscores %>% select(Chromosome, Position, GenecodeComprehensiveInfo, "CaddPhred","ApcProteinFunctionV3", "ApcConservationV2", "ApcEpigeneticsActive", "ApcEpigeneticsRepressed",
                                 "ApcEpigeneticsTranscription", "ApcLocalNucleotideDiversityV3", "ApcMutationDensity",
                                 "ApcTranscriptionFactor", "ApcMappability")
fwrite(inscores.aPCs_df,paste0("output/fnc_annot_inscores",region),sep="\t",quote=F,row.names=F)

#Clinical annotation:
clinvar <- favor.digest %>% filter(!is.na(Clnsig))
clin_genes <- unique(clinvar$Geneinfo)
clinvar_df <- clinvar %>% select(Chromosome, Position, Clnsig, Geneinfo)
fwrite(clinvar_df,paste0("output/fnc_annot_clinvar",region),sep="\t",quote=F,row.names=F)

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
fwrite(as.data.frame(nearest_gene),paste0("input/nearest_genes_raw",region),row.names=FALSE, col.names=FALSE,quote=F)
fwrite(varannot_genes,paste0("/alice-home/3/n/nnp5/PhD/PhD_project/Var_to_Gene/input/varannot_genes_",region),row.names=FALSE, col.names=FALSE,quote=F)
#create /alice-home/3/n/nnp5/PhD/PhD_project/Var_to_Gene/input/var2genes_raw",region,".xlsx" and add nearest_genes and varannot_genes sheets.