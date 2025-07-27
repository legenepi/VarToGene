#!/usr/bin/env Rscript

#Rationale: visualisation of genes prioritised from all analyses with updated credset and V2G results as March 2025.


library(scales)
library(tidyverse)
library(data.table)
library(readxl)
library(RColorBrewer)


#input each analysis table:
#genes_raw <- read_table("output/V2G_March2025/Loci2Gene_March2025.txt")
genes_raw <- read_table("output/V2G_March2025/Loci2Gene_May2025_nearestgeneOTP.txt")
v2g_full <- unique(expand.grid(x = genes_raw$Gene, y = genes_raw$Evidence, KEEP.OUT.ATTRS = TRUE)) %>%
            arrange(x) %>% rename("Gene" = x, "Evidence" = y)
mask <- as.data.frame(do.call(paste0, v2g_full) %in% do.call(paste0, genes_raw[,-1])) %>% rename(status='do.call(paste0, v2g_full) %in% do.call(paste0, genes_raw[, -1])')
v2g_full_combination <- cbind(v2g_full,mask) %>% mutate(status=as.integer(as.logical(status))) %>% arrange(status,Gene)

locus_v2g_full_combination <- v2g_full_combination %>% left_join(genes_raw, by = "Gene", relationship = "many-to-many") %>%
                              select(-Evidence.y) %>% rename(Evidence = "Evidence.x") %>% unique()

#Do heatmap for each Locus:
#fnct for the plot:
fc_heatmap_all <- function(df,x_val,y_val,fill_val,locus_par) {
      ggplot(df, aes(x=x_val, y=y_val, fill=fill_val)) +
      geom_tile() +
      scale_fill_gradient2(low = "white",
                       mid = "white",
                       high = "#9EBCDA") +
      scale_x_discrete(position = "top") +
      theme(axis.text = element_text(face = "bold"), axis.ticks.y=element_blank(), axis.ticks.x=element_blank()) +
      xlab("Evidence") + ylab("Gene") + theme(legend.position = "none") +
      ggtitle(locus_par)
      }

loci <- locus_v2g_full_combination %>% select(Locus) %>% unique()

for (i in seq(1:dim(loci)[1])) {
#locus_var <- "10_rs201499805_rs1444789_8542744_9564361"
locus_var <- loci[i,]
#locus <- as.character(locus_v2g_full_combination %>% filter(Locus == locus_var) %>% select(Locus) %>% unique())
df <- locus_v2g_full_combination %>% filter(Locus == locus_var) %>% select(-Locus) %>% arrange(Gene, Evidence)
#png(paste0("./output/Additional_credset_snps_March2025_output/V2G_heatmap_",as.character(locus_var),".png"), width=1000, height = 950, res = 150)
plot_df <- fc_heatmap_all(df,df$Evidence,df$Gene,df$status,locus_var) +
        theme(strip.background = element_blank(), strip.text = element_blank(),
        axis.text.x=element_text(angle=75, vjust=0, hjust=0, face="bold"),
        axis.title=element_blank(),
        axis.text.y=element_text(hjust=1, face="italic"),
        axis.ticks = element_blank(),
        legend.title = element_blank()) +
        coord_fixed(ratio = 1) + coord_equal()
#print(plot_df)
#dev.off()
ggsave(plot_df, file=paste0("./output/Additional_credset_snps_March2025_output/V2G_heatmap_",as.character(locus_var),".png"),  dpi=90, width = 5, height = 7)
}
