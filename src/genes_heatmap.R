#!/usr/bin/env Rscript

#Rationale: visualisation of genes prioritised from all analyses

library(scales)
library(tidyverse)
library(data.table)
library(readxl)
library(RColorBrewer)


#input each analysis table:
genes_raw <- "input/var2genes_raw.xlsx"
nearest_gene <- read_excel(genes_raw, sheet = "nearest_genes",col_names = "gene")
nearest_gene$evidence <- "nearest"
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
fwrite(v2g_full_combination2,"./output/v2g_gene_prioritisation.txt",sep="\t",quote=F)

# Plot
##use this to find which colour I want: fro the pie, I extract the index for the colour I want:
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


#Split plot into 9 plots of 11 genes each:
test <- v2g_full_combination2 %>%
    arrange(gene, variable) %>%
    group_by() %>%
    mutate(facet=c(rep(9, ceiling(n()/9)),rep(8, ceiling(n()/9)),rep(7, ceiling(n()/9)),rep(6, ceiling(n()/9)),rep(5, ceiling(n()/9)),rep(4, ceiling(n()/9)),rep(3, ceiling(n()/9)),rep(2, ceiling(n()/9)), rep(1, floor(n()/9)))) %>%
    ungroup

png("./output/V2G_heatmap_subplots.png", width=1000, height = 900)
fc_heatmap_all(test,test$variable,test$gene,test$value) + facet_wrap(~facet, scales="free", ncol=3) +
        theme(strip.background = element_blank(), strip.text = element_blank(),
        axis.text.x=element_text(angle=75, vjust=0, hjust=0, face="bold"),
        axis.title=element_blank(),
        axis.text.y=element_text(hjust=1, face="italic"),
        axis.ticks = element_blank(),
        legend.title = element_blank())
dev.off()