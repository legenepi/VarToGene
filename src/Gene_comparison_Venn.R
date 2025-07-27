#!/usr/bin/env Rscript

#Rationale: Venn Diagram comparison FINEMAP-SuSie VS SuSie finemapping and variant-to-gene mapping.

# Load library
library(tidyverse)
library(data.table)
library(VennDiagram)

# input data
genes <- read_table("input/Additional_credset_snps_March2025/Genes_comparison.txt")
finemapsusie <- genes$Finemap_SuSiE
susie <- genes$SuSie



# Chart
venn.diagram(
  x = list(finemapsusie, susie),
  category.names = c("finemapsusie", "susie"),
  filename = 'output/Additional_credset_snps_March2025_output/Genes_comparison_Venn.png',
  output=TRUE
)