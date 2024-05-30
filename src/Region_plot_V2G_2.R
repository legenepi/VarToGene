#!/usr/bin/env Rscript


#Rationale: plot of fine-mapping PIP with functional annotation for variants in credible set:

library(tidyverse)
library(data.table)
library(dplyr)
library(RColorBrewer)
library(cowplot)
library(readxl)

args <- commandArgs(T)

r2_file <- args[1]
locus <- as.character(args[2])
locus_ggtitle <- as.character(args[3])
start <- as.numeric(args[4])
end <- as.numeric(args[5])
chr <- args[6]
snp_lead <- as.character(args[7])
print(start)
print(end)
print(chr)
print(snp_lead)
#cp /scratch/gen1/nnp5/Var_to_Gen_tmp/regionalplot/gwas_credset_annot input/gwas_credset_annot
#cp /scratch/gen1/nnp5/Var_to_Gen_tmp/regionalplot/gwas_credset_annot_chr3_49524027_50524027_rs778801698 input/gwas_credset_annot_chr3_49524027_50524027_rs778801698
df <- fread("input/gwas_credset_annot") %>% select(-snpid)
#if chr3 s778801698:
df <- fread("input/gwas_credset_annot_chr3_49524027_50524027_rs778801698") %>% select(-snpid)
df <- df %>% mutate(snpid = paste0(chromosome,"_",posb37,"_",pmin(a1,a2),"_",pmax(a1,a2)))
df_sugg <- df %>% filter(pval <= 0.000005)
r2 <- fread(r2_file) %>% select(CHR_B,BP_B,SNP_B,R2)
colnames(r2) <- c("chromosome","posb37","snpid","R2")

df_r2 <- left_join(df,r2,by=c("chromosome","posb37","snpid"))

df_r2 <- df_r2 %>% mutate(Functional_annotation=ifelse(is.na(df_r2$Functional_annotation), "unknown",df_r2$Functional_annotation))
df_r2$Functional_annotation <- as.factor(df_r2$Functional_annotation)
df_r2$logP <- -log10(df$pval)

#for R2 color gradient in the plot:
myPalette <- colorRampPalette(rev(brewer.pal(11, "Spectral")))
sc <- scale_colour_gradientn(colours = myPalette(100), limits=c(0, 1))

##Plot with functional annotation and R2:
annot_r2 <- df_r2 %>% ggplot(aes(posb37, logP, shape = Functional_annotation, colour = R2, label = snpid)) +
  geom_point(size=2) +
  scale_shape_manual(values=c("downstream"=15, "exonic"=16, "intergenic"=17, "intronic"=18, "ncRNA_exonic"=19, "ncRNA_intronic"=10, "unknown"=11,"upstream;downstream"=12, "UTR3"=13, "UTR5"=14)) +
  geom_point(data=df_r2[df_r2$credset==1,], pch=21, fill=NA, size=4, colour="red", stroke=1, alpha=0.5) +
  #geom_text(aes(label=ifelse(R2>0.95,as.character(snpid),'')),hjust=0,vjust=0,size=2.5,colour="black") + #the rsid names create messy plot - removed it
  theme_minimal() +
  theme(legend.position = "right") + sc +
  ggtitle(paste0(locus_ggtitle," (snp lead:", snp_lead,")")) + theme(plot.title = element_text(hjust=0.5)) + xlim(start,end) + ylab("-log10(pvalue)")


##Plot with functional annotation, R2 and gene location:
#Treated the gene plot as a gantt chart plot !
#prioritised genes and evidence:
#l2g <- read_excel("src/report/locus2gene.xlsx",sheet = "L2G_clean") %>% filter(grepl(end,locus)) %>% select(gene)
#if chr3 rs778801698:
l2g <- read_excel("src/report/locus2gene.xlsx",sheet = "L2G_clean_chr3_rs778801698") %>% filter(grepl(end,locus)) %>% select(gene) %>% unique()

#This gene list is from the R package LocusZooms
genes <- fread("/home/n/nnp5/software/LocusZooms/Gencode_GRCh37_Genes_UniqueList2021.txt")
#filter gene in the locus:
genes <- genes %>% filter(Chrom==paste0("chr",chr), Start >= start, End <= end)

## Set factor level to order the genes on the plot
genes$Gene <- as.factor(genes$Gene)
plot_gantt_gene <- qplot(ymin = Start,
                    ymax = End,
                    x = Gene,
                    colour = Coding,
                    geom = "linerange",
                    data = genes,
                    size = I(5)) +
    scale_colour_manual(values = c("lncRNA"="burlywood3", "ncRNA"="aquamarine4", "proteincoding"="chocolate", "pseudogene"="chocolate4")) +
    coord_flip() +
    theme_bw() +
    theme(panel.grid = element_blank()) +
    ylab("posb37") +
    xlab("genes") +
    ylim(start,end) + labs(title = paste0("V2G identified genes: ",l2g))


plot_grid(annot_r2 + theme(legend.justification = c(0,1)), plot_gantt_gene + theme(legend.justification = c(0,1)), ncol=1, align='v', rel_heights = c(2,1))
ggsave(locus_ggtitle, width = 80, height = 60, units = "cm")
ggsave(locus_ggtitle, width = 50, height = 60, units = "cm")