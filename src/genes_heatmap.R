#!/usr/bin/env Rscript

#Rationale: visualisation of genes prioritised from all analyses
#NB: add genes related to chr 3 rs778801698 and delete genes of the MHC region.

library(scales)
library(tidyverse)
library(data.table)
library(readxl)
library(RColorBrewer)


#input each analysis table:
genes_raw <- "input/var2genes_raw.xlsx"
nearest_gene <- read_excel(genes_raw, sheet = "nearest_genes",col_names = "gene")
nearest_gene$evidence <- "nearest"
#Edit for locus 5_rs2188962_rs848_131270805_132496500 with sentinel 5_131885240_G_C: checked on NCBI, it is an intronic
#variant of IL5, so nearest gene is IL5; remove AC116366.3:
nearest_gene <- nearest_gene %>% filter(gene != "AC116366.3")
annotation <- read_excel(genes_raw, sheet = "varannot_genes",col_names = "gene")
annotation$evidence <- as.factor("func_annot")
eqtl <- read_excel(genes_raw, sheet = "eQTL_genes_merge")
pqtl <- read_excel(genes_raw, sheet = "pQTL_genes_merge")
pops <- read_excel(genes_raw, sheet = "PoPS_genes",col_names = "gene")
pops$evidence <- as.factor("PoPS")
mouseko <- read_excel(genes_raw, sheet = "Mouseko_genes",col_names = "gene")
mouseko$evidence <- as.factor("mouse_KO")
raredis <- read_excel(genes_raw, sheet = "Raredisease_genes",col_names = "gene")
raredis$evidence <- as.factor("rare_disease")
#merge:
v2g_full <- rbind(nearest_gene,annotation,eqtl,pqtl,pops,mouseko,raredis)
#Rename ILRL1 and ST2 into unique name: "IL1RL1/ST2"
v2g_full$gene <- sub("^IL1RL1$","IL1RL1/ST2",v2g_full$gene)
v2g_full$gene <- sub("^ST2$","IL1RL1/ST2",v2g_full$gene)
#genes with mhc and without chr3:
length(unique(v2g_full$gene))

#chr3: do not indlude anymore as it did not replicate without SAPR I-II and SARP III
#genes_raw_chr3 <- "input/var2genes_raw_chr3_49524027_50524027_rs77880169.xlsx"
#nearest_gene_chr3 <- read_excel(genes_raw_chr3, sheet = "nearest_genes",col_names = "gene")
#nearest_gene_chr3$evidence <- "nearest"
#annotation_chr3 <- read_excel(genes_raw_chr3 , sheet = "varannot_genes",col_names = "gene")
#annotation_chr3$evidence <- as.factor("func_annot")
#eqtl_chr3 <- read_excel(genes_raw_chr3, sheet = "eQTL_genes_merge",col_names = "gene")
#eqtl_chr3$evidence <- "eQTL"
#pqtl_chr3 <- read_excel(genes_raw_chr3, sheet = "pQTL_genes_merge",col_names = "gene")
#pqtl_chr3$evidence <- "pQTL"
#pops_chr3 <- read_excel(genes_raw_chr3, sheet = "PoPS_genes",col_names = "gene")
#pops_chr3$evidence <- as.factor("PoPS")
#mouseko_chr3 <- read_excel(genes_raw_chr3, sheet = "Mouseko_genes",col_names = "gene")
#mouseko_chr3$evidence <- as.factor("mouse_KO")
#raredis_chr3 <- read_excel(genes_raw_chr3, sheet = "Raredisease_genes",col_names = "gene")
#raredis_chr3$evidence <- as.factor("rare_disease")
#merge:
#v2g_full_chr3 <- rbind(nearest_gene_chr3,annotation_chr3,eqtl_chr3,pqtl_chr3,pops_chr3,mouseko_chr3,raredis_chr3)
#genes chr3:
#length(unique(v2g_full_chr3$gene))

#merge previous results with chr3 results:
#v2g_full <- rbind(v2g_full, v2g_full_chr3)
#genes with mhc and with chr3:
#length(unique(v2g_full$gene))

#delete genes of MHC region:
mhc_locus_gene <- read_excel("/home/n/nnp5/PhD/PhD_project/Var_to_Gene/src/report/var2gene_full.xlsx", sheet = "Sheet1") %>%
                  filter(locus == "6_rs9271365_32086794_33086794")
mhc_genes <- unique(gsub("\\(.*", "",mhc_locus_gene$gene))
v2g_full <- v2g_full %>% filter(! gene %in% mhc_genes)
#genes without mhc and without chr3:
length(unique(v2g_full$gene))

#table with all the possible combination:
v2g_full_combination <- unique(expand.grid(x = v2g_full$gene, y = v2g_full$evidence, KEEP.OUT.ATTRS = TRUE)) %>% arrange(x)
colnames(v2g_full_combination) <- colnames(v2g_full)
#mask (conceptually: v2g_full_combination %in% v2g_full):
mask <- as.data.frame(do.call(paste0, v2g_full_combination) %in% do.call(paste0, v2g_full)) %>% rename(status='do.call(paste0, v2g_full_combination) %in% do.call(paste0, v2g_full)')
v2g_full_combination <- cbind(v2g_full_combination,mask) %>% mutate(status=as.integer(as.logical(status))) %>% arrange(status,gene)
#long to wide reshape using spread() in tidyr:
v2g_full_combination <- spread(v2g_full_combination, evidence, status) %>% remove_rownames %>% column_to_rownames(var="gene")

#  Convert row names into first column
v2g_full_combination <- setDT(v2g_full_combination, keep.rownames = TRUE)[]
setnames(v2g_full_combination, "rn", "gene")

# reshape your data
v2g_full_combination2 <- melt(v2g_full_combination, id.var = "gene")
length(unique(v2g_full_combination2$gene))
#fwrite(v2g_full_combination2,"./output/v2g_gene_prioritisation_chr3_noMHC.txt",sep="\t",quote=F)
#v2g_full_combination2 <- fread("./output/v2g_gene_prioritisation_chr3_noMHC.txt")
fwrite(v2g_full_combination2,"./output/v2g_gene_prioritisation_nochr3_noMHC.txt",sep="\t",quote=F)
v2g_full_combination2 <- fread("./output/v2g_gene_prioritisation_nochr3_noMHC.txt")



# Plot
##use this to find which colour I want: for the pie, I extract the index for the colour I want:
##everytime it changes, so I noted down the actual code of the colour: ##18: col[18] #9EBCDA
#n <- 30
#colrs <- brewer.pal.info[brewer.pal.info$colorblind == TRUE, ]
#col_vec = unlist(mapply(brewer.pal, colrs$maxcolors, rownames(colrs)))
#col <- sample(col_vec, n)
#area <- rep(1,n)
#pie(area, col = col)

#fnct for the plot:
fc_heatmap_all <- function(df,x_val,y_val,fill_val) {
      ggplot(df, aes(x=x_val, y=y_val, fill=fill_val)) +
      geom_tile() +
      scale_fill_gradient2(low = "white",
                       mid = "white",
                       high = "#9EBCDA") +
      scale_x_discrete(position = "top") +
      theme(axis.text = element_text(face = "bold"), axis.ticks.y=element_blank(), axis.ticks.x=element_blank()) +
      xlab("Evidence") + ylab("Gene") + theme(legend.position = "none") +
      ggtitle("Gene prioritization")
      }


#Split plot into 3 plots (tot of 89 genes) HAVING PROBLEMS, NEED TO RESOLVE
#data handling - add an empty row to have 90 rows and split the heatmap by 3 plots of 30 genes (rows):
df<-data.frame(c("  ","  ","  ","  ","  ","  ","  "),as.factor(c("nearest","func_annot","eQTL","pQTL","PoPS","mouse_KO","rare_disease")),c(0,0,0,0,0,0,0))
names(df)<-c("gene","variable","value")
df$value <- as.numeric(df$value)
v2g_full_combination3 <- rbind(v2g_full_combination2,df)
v2g_full_combination3$variable <- ordered(v2g_full_combination3$variable, levels = c("nearest","func_annot","eQTL","pQTL","PoPS","mouse_KO","rare_disease"))

test <- v2g_full_combination3 %>%
    arrange(gene, variable) %>%
    group_by() %>%
    mutate(facet=c(rep(3, ceiling(n()/3)),rep(2, ceiling(n()/3)),rep(1, ceiling(n()/3)))) %>%
    ungroup

png("./output/V2G_heatmap_subplots_nochr3_noMHC.png", width=1000, height = 950)
fc_heatmap_all(test,test$variable,test$gene,test$value) + facet_wrap(~facet, scales="free", ncol=3) +
        theme(strip.background = element_blank(), strip.text = element_blank(),
        axis.text.x=element_text(angle=75, vjust=0, hjust=0, face="bold"),
        axis.title=element_blank(),
        axis.text.y=element_text(hjust=1, face="italic"),
        axis.ticks = element_blank(),
        legend.title = element_blank())
dev.off()


#Split plot into 3 plots of 37 genes each (tot of 111 genes): this was with chr3 genes - deprecated
test <- v2g_full_combination2 %>%
    arrange(gene, variable) %>%
    group_by() %>%
    mutate(facet=c(rep(3, ceiling(n()/3)),rep(2, ceiling(n()/3)),rep(1, ceiling(n()/3)))) %>%
    ungroup

png("./output/V2G_heatmap_subplots_chr3_noMHC.png", width=1000, height = 950)
fc_heatmap_all(test,test$variable,test$gene,test$value) + facet_wrap(~facet, scales="free", ncol=3) +
        theme(strip.background = element_blank(), strip.text = element_blank(),
        axis.text.x=element_text(angle=75, vjust=0, hjust=0, face="bold"),
        axis.title=element_blank(),
        axis.text.y=element_text(hjust=1, face="italic"),
        axis.ticks = element_blank(),
        legend.title = element_blank())
dev.off()


#Additional:
#Edit heatmap so that the genes are clustered by genomic loci:
#fill_val needs to be a color for each locus.
#Create a custom color scale
library(RColorBrewer)
myColors <- brewer.pal(17,"Set1")
names(myColors) <- levels(test$locus)
colScale <- scale_colour_manual(name = "grp",values = myColors)

#Split plot into 3 plots of 37 genes each (tot of 111 genes):
gene_locus <- read_excel("src/report/locus2gene.xlsx", sheet = "L2G_clean_noMHC_chr3")
v2g_locus <- v2g_full_combination2 %>% left_join(gene_locus, by = "gene")
v2g_locus <- value

test <- v2g_full_combination2 %>%
    arrange(gene, variable) %>%
    group_by() %>%
    mutate(facet=c(rep(3, ceiling(n()/3)),rep(2, ceiling(n()/3)),rep(1, ceiling(n()/3)))) %>%
    ungroup

test2 <- v2g_locus %>%
    arrange(locus, gene, variable) %>%
    group_by() %>%
    mutate(facet=c(rep(3, ceiling(n()/3)),rep(2, ceiling(n()/3)),rep(1, ceiling(n()/3)))) %>%
    ungroup


locus_unique <- as.data.frame(levels(as.factor(test2$locus)))
colnames(locus_unique) <- "locus"
locus_unique$index <- seq(1,17)
test2 <- test2 %>% left_join(locus_unique, by = "locus")
test2 <- test2 %>% mutate(value_locus = ifelse(value == 1, index, 0))
myColors <- c(brewer.pal(8,"Pastel2"),brewer.pal(9,"Pastel1"))

fc_heatmap_all <- function(df,x_val,y_val,fill_val) {
      ggplot(df, aes(x=x_val, y=y_val, fill=fill_val)) +
      geom_tile() +
      scale_fill_manual(breaks = levels(fill_val), values = myColors) +
      scale_x_discrete(position = "top") +
      theme(axis.text = element_text(face = "bold"), axis.ticks.y=element_blank(), axis.ticks.x=element_blank()) +
      xlab("Evidence") + ylab("Gene") + theme(legend.position = "none") +
      ggtitle("Gene prioritization")
      }


plot <- fc_heatmap_all(test2,test2$variable,test2$gene,test2$value_locus) + facet_wrap(~facet, scales="free", ncol=3) +
        theme(strip.background = element_blank(), strip.text = element_blank(),
        axis.text.x=element_text(angle=75, vjust=0, hjust=0, face="bold"),
        axis.title=element_blank(),
        axis.text.y=element_text(hjust=1, face="italic"),
        axis.ticks = element_blank(),
        legend.title = element_blank())


plot
