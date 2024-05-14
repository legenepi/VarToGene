#!/usr/bin/env Rscript

#Rationale: Variant Annotation with FAVOR - need a new format because the colnames have changed as FAVOR got updated

sink(stderr())

suppressMessages(library(tidyverse))
suppressMessages(library(data.table))
suppressMessages(library(corrplot))
#suppressMessages(library(xlsx))
#suppressMessages(library(qdap))
suppressMessages(library(VennDiagram))
library(RColorBrewer)

#args <- commandArgs(T)
favor_file <- args[1]
region <- args[2]
#favor_file <- "/alice-home/3/n/nnp5/PhD/PhD_project/Var_to_Gene/input/FAVOR_credset_chrpos38_2024_05_14_chr3_49524027_50524027_rs778801698.txt.csv"
favor <- fread(favor_file)

favor.digest <- favor %>% select("VariantVcf", "Chromosome", "Position", "Rsid", "GenecodeComprehensiveCategory",
                                 "GenecodeComprehensiveInfo", "GenecodeComprehensiveExonicCategory", "GenecodeComprehensiveExonicInfo",
                                 "CAGEPromoter", "CAGEEnhancer", "Genehancer", "SuperEnhancer", "CaddPhred", "FathmmXf",
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
locus_sig_list <- sig_list_tmp %>% filter(locus == as.character(i))
locus_sig_list <- locus_sig_list %>% filter(PIP_average == max(locus_sig_list$PIP_average)) %>% select(locus,sentinel,chr,posb37)
sentinels <- rbind(sentinels,locus_sig_list)
seninels <- sentinels %>% filter(locus == "3_rs778801698_49524027_50524027")

#liftover fo rhte senitnel to find b38:
posb38 <- fread("input/highestPIP_sentinels_hglft_genome_b38.bed") %>% rename(posb38=V2)
posb38  <- posb38 %>% separate(V4, c("chr", "pos2", "posb37"))
posb38  <- posb38 %>% select(chr, posb37, posb38)
posb38$chr <- as.numeric(gsub("chr","",posb38$chr))
posb38$posb37 <- as.numeric(posb38$posb37)
sentinel_b38 <- sentinels %>% left_join(posb38, by =c("chr", "posb37")) %>% rename(Chromosome=chr, Position=posb38)



###Nearest genes:
favor.digest.NG <- favor.digest %>% select(Chromosome, Position, "Genecode.Comprehensive.Info") %>% rename(Nearest_gene="Genecode.Comprehensive.Info")
sentinel_b38_NG <- sentinel_b38 %>% left_join(favor.digest.NG, by = c("Chromosome", "Position")) %>% rename(posb38=Position)
#Take closet genes when multiples are annotated:
#1:AP001189.5
#5:TSLP
#8:IL33
#9:IL1RL1
#16:HLA-DQA1
#17:AC044784.1
sentinel_b38_NG$Nearest_gene[1] <- "AP001189.5"
sentinel_b38_NG$Nearest_gene[5] <- "TSLP"
sentinel_b38_NG$Nearest_gene[8] <- "IL33"
sentinel_b38_NG$Nearest_gene[9] <- "IL1RL1"
sentinel_b38_NG$Nearest_gene[16] <- "HLA-DQA1"
sentinel_b38_NG$Nearest_gene[17] <- "AC044784.1"
fwrite(sentinel_b38_NG,paste0("output/PIP_sentinels_nearestgenes",region),sep="\t",quote=F)
nearest_gene <- unique(unlist(strsplit(sentinel_b38_NG$Nearest_gene,",")))

#Categorical annotation:
##Exploratory:
summary(as.factor(favor.digest$"Genecode.Comprehensive.Category"))
summary(as.factor(favor.digest$"Genecode.Comprehensive.Info"))
##Info for the 5 coding variants:
summary(as.factor(favor.digest$"Genecode.Comprehensive.Exonic.Category"))
summary(as.factor(favor.digest$"Genecode.Comprehensive.Exonic.Info"))

#FAVOR functional:
#FANTOM5 CAGE promoter or enhancer, Functional Integrative score > 15, ClinVar
#CAGE promoter or enhancer (from FANTOM5):
fantom5 <- favor.digest %>% filter(!is.na(CAGE.Promoter) | !is.na(CAGE.Enhancer))
fantom5_genes <- unique(fantom5$Genecode.Comprehensive.Info)
fantom5_df <- fantom5 %>% select(Chromosome, Position, "Genecode.Comprehensive.Info", CAGE.Promoter, CAGE.Enhancer) %>% rename(Nearest_gene="Genecode.Comprehensive.Info")
fwrite(fantom5_df,paste0("output/fnc_annot_fantom5",sep="\t",region)quote=F,row.names=F)


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
#Save the corplot x the report:
png(file = "/home/n/nnp5/PhD/PhD_project/Var_to_Gene/output/corrplot_integrative_score_15.png")
corplot_integrativescore <- corrplot(M2, method="number")
dev.off()

inscores.aPCs_df <- inscores %>% select(Chromosome, Position, Genecode.Comprehensive.Info, "CADD.phred","aPC.Protein.Function", "aPC.Conservation","aPC.Epigenetics.Active",
                                      "aPC.Epigenetics.Repressed","aPC.Epigenetics.Transcription","aPC.Local.Nucleotide.Diversity",
                                      "aPC.Mutation.Density", "aPC.Transcription.Factor", "aPC.Mappability")
fwrite(inscores.aPCs_df,paste0("output/fnc_annot_inscores",sep="\t",region),quote=F,row.names=F)

#Clinical annotation:
clinvar <- favor.digest %>% filter(!is.na(Clinical.Significance))
clin_genes <- unique(clinvar$Gene.Reported)
clinvar_df <- clinvar %>% select(Chromosome, Position, Clinical.Significance, Gene.Reported)
fwrite(clinvar_df,paste0("output/fnc_annot_clinvar",sep="\t",region),quote=F,row.names=F)

#Genes from Variant Annotation:
#fantom5_genes, inscores_genes, clin_genes
##Polish the genes (remove GeneID, remove ensembl, remove where mutation occurs, keep only gene name)
#fantom_genes: everything inside parenthesis needs to be deleted, and divide genes listed together
fantom5_genes <- bracketX(fantom5_genes) %>% unique()
fantom5_genes <- unlist(strsplit(fantom5_genes, "\\s*,\\s*"))

#inscores_genes: everything inside parenthesis needs to be deleted, and divide genes listed together
inscores_genes <- bracketX(inscores_genes) %>% unique()
inscores_genes <- unlist(strsplit(inscores_genes, "\\s*,\\s*"))
inscores_genes <- unlist(strsplit(inscores_genes, ';', fixed=TRUE))

#clin_genes: separate element if '|' present, and then remove GeneID-after the ':'
clin_genes <- unlist(strsplit(clin_genes, '|', fixed=TRUE))
clin_genes <- gsub(":.*", "", clin_genes) %>% unique()

#Save the genes:
fwrite(as.data.frame(nearest_gene),paste0("input/nearest_genes_raw",region),row.names=FALSE, col.names=FALSE,quote=F)
write.xlsx(varannot_genes,paste0("/alice-home/3/n/nnp5/PhD/PhD_project/Var_to_Gene/input/var2genes_raw",region,".xlsx"),sheetName = "varannot_genes", row.names=FALSE, col.names=FALSE)