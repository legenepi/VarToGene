#!/usr/bin/env Rscript

#Rationale: Variant Annotation with FAVOR

#Run as:
#Rscript src/ \


sink(stderr())

suppressMessages(library(tidyverse))
suppressMessages(library(data.table))
suppressMessages(library(corrplot))
suppressMessages(library(xlsx))

args <- commandArgs(T)

favor_file <- args[1]
favor_file <- "/alice-home/3/n/nnp5/PhD/PhD_project/Var_to_Gene/input/FAVOR_credset_chrpos38_2023_08_08.csv"
favor <- fread(favor_file)

favor.digest <- favor %>% select("Variant (VCF)", "Chromosome", "Position", "rsID", "Genecode Comprehensive Category",
                                 "Genecode Comprehensive Info", "Genecode Comprehensive Exonic Category", "Genecode Comprehensive Exonic Info",
                                 "CAGE Promoter", "CAGE Enhancer", "GeneHancer", "SuperEnhancer", "CADD phred", "Fathmm XF",
                                 "aPC-Protein-Function", "aPC-Conservation", "aPC-Epigenetics-Active", "aPC-Epigenetics-Repressed",
                                 "aPC-Epigenetics-Transcription", "aPC-Local-Nucleotide-Diversity", "aPC-Mutation-Density",
                                 "aPC-Transcription-Factor", "aPC-Mappability", "Clinical Significance",
                                 "Clinical Significance (genotype includes)", "Disease Name",
                                 "Disease Name (included variant)", "Review Status", "Allele Origin",
                                 "Disease Database ID", "Disease Database ID (included variant)", "Gene Reported")

colnames(favor.digest) <- make.names(colnames(favor.digest), unique=TRUE)

#Categorical annotation:
##Exploratory:
summary(as.factor(favor.digest$"Genecode.Comprehensive.Category"))
summary(as.factor(favor.digest$"Genecode.Comprehensive.Info"))
##Info for the 5 coding variants:
summary(as.factor(favor.digest$"Genecode.Comprehensive.Exonic.Category"))
summary(as.factor(favor.digest$"Genecode.Comprehensive.Exonic.Info"))

#CAGE promoter or enhancer (from FANTOM5):
fantom5 <- favor.digest %>% filter(!is.na(CAGE.Promoter) | !is.na(CAGE.Enhancer))
fantom5_genes <- unique(fantom5$Genecode.Comprehensive.Info)


#Functional integrative score:
#Among 562 (53 removed because NAs) there are a correlations between aPC.Epigenetics.Active and aPC.Transcription.Factor (0.75), between aPC.Mutation.Density
#and aPC.Local.Nucleotide.Diversity (0.80); between CADD.phred and aPCs.Conservation (0.83).
favor.aPCs <- favor.digest %>% select("CADD.phred","aPC.Protein.Function", "aPC.Conservation","aPC.Epigenetics.Active",
                                      "aPC.Epigenetics.Repressed","aPC.Epigenetics.Transcription","aPC.Local.Nucleotide.Diversity",
                                      "aPC.Mutation.Density", "aPC.Transcription.Factor", "aPC.Mappability")
favor.aPCs <- na.omit(favor.aPCs)
M <- cor(favor.aPCs)
head(round(M,3))
corrplot(M, method="number")

##Integrative score > 15:
inscores <- favor.digest %>% filter(CADD.phred > 15 | aPC.Protein.Function > 15 | aPC.Conservation > 15 | aPC.Epigenetics.Active > 15 | aPC.Epigenetics.Repressed > 15 |
aPC.Epigenetics.Transcription > 15 | aPC.Local.Nucleotide.Diversity > 15 | aPC.Mutation.Density > 15 | aPC.Transcription.Factor > 15 | aPC.Mappability > 15)
inscores_genes <- unique(inscores$Genecode.Comprehensive.Info)

#Among 99 integrative score > 15 variants, these correlations change a little: between aPC.Epigenetics.Active and aPC.Transcription.Factor (0.69), between aPC.Mutation.Density
#and aPC.Local.Nucleotide.Diversity (0.63); between CADD.phred and aPCs.Conservation (0.91).
inscores.aPCs <- inscores %>% select("CADD.phred","aPC.Protein.Function", "aPC.Conservation","aPC.Epigenetics.Active",
                                      "aPC.Epigenetics.Repressed","aPC.Epigenetics.Transcription","aPC.Local.Nucleotide.Diversity",
                                      "aPC.Mutation.Density", "aPC.Transcription.Factor", "aPC.Mappability")
M2 <- cor(inscores.aPCs)
head(round(M,3))
corrplot(M2, method="number")

#Clinical annotation:
clinvar <- favor.digest %>% filter(!is.na(Clinical.Significance))
clin_genes <- unique(clinvar$Gene.Reported)

#Genes from Variant Annotation:
#fantom5_genes, inscores_genes, clin_genes
varannot_genes <- list(unique(c(fantom5_genes,inscores_genes,clin_genes)))

#Save the genes:
write.xlsx(varannot_genes,"/alice-home/3/n/nnp5/PhD/PhD_project/Var_to_Gene/input/var2genes_raw.xlsx",sheetName = "varannot_genes", row.names=FALSE)

