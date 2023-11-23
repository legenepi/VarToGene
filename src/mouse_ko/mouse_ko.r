#!/usr/bin/env Rscript

#Rationale: Find nearby human ortholog of mouse knockout genes that are within +/- 500Kb from SA GWAS sentinels.

library(GenomicRanges)
library(rtracklayer)
library(tidyverse)
library(pander)
library(scales)
library(readxl)
library(data.table)

setwd("/scratch/gen1/nnp5/Var_to_Gen_tmp/mouse_ko/")

DISTANCE <- 5e5

getN <- function(x, y) x %>% pull({{y}}) %>% n_distinct


##Retrieve Sentinel SNPs - highest PIP in the fine-mapped loci:
sig_list_tmp <- fread("/data/gen1/UKBiobank_500K/severe_asthma/Noemi_PhD/data/replsugg_valid_credset.txt")
setnames(sig_list_tmp,"chromosome","chr")
setnames(sig_list_tmp,"position","pos")
sig_list_tmp$sentinel <- paste0(sig_list_tmp$chr,"_",sig_list_tmp$pos,"_",sig_list_tmp$allele1,"_",sig_list_tmp$allele2)
locus <- unique(sig_list_tmp$locus)

sentinels <- data.frame(matrix(ncol = 4,nrow = 0))
colnames(sentinels) <- c("locus","sentinel","chr","pos")
for(i in locus){
    locus_sig_list <- sig_list_tmp %>% filter(locus == as.character(i))
    locus_sig_list <- locus_sig_list %>% filter(PIP_average == max(locus_sig_list$PIP_average)) %>% select(locus,sentinel,chr,pos)
    sentinels <- rbind(sentinels,locus_sig_list)
    }

## retrieve genomic reference for genes
#resolved problem:
#if (!require("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install(version = "3.18")

ucsc <- browserSession("UCSC")
genome(ucsc) <- "hg19"
refseq <- ucsc %>%
  ucscTableQuery(track="NCBI RefSeq", table="refGene") %>%
  getTable %>%
  as_tibble

## Mouse genotype-phenotype data
## Store in /scratch/gen1/nnp5/Var_to_Gen_tmp/mouse_ko/ and downloaded from:
#latest release 2023-07-06:
#wget -P ${tmp_path}/mouse_ko/ https://ftp.ebi.ac.uk/pub/databases/impc/all-data-releases/latest/results/genotype-phenotype-assertions-ALL.csv.gz
#latest release 2023-11-22:
#wget -P ${tmp_path}/mouse_ko/ http://ftp.ebi.ac.uk/pub/databases/genenames/hcop/human_mouse_hcop_fifteen_column.txt.gz

orthologs <- read_tsv("human_mouse_hcop_fifteen_column.txt.gz")
all_geno_pheno <- read_csv("genotype-phenotype-assertions-ALL.csv.gz")

geno_pheno <- all_geno_pheno %>%
  separate(top_level_mp_term_name, c("top_level_mp_term", "top_level_mp_subtype"), ",")

#display all the possible top level mp term and related number of genes:
geno_pheno %>%
  group_by(top_level_mp_term) %>%
  summarise(nlines=n(), ngenes=n_distinct(marker_symbol)) %>%
  drop_na %>%
  pander(big.mark=",")


## Overlap of severe asthma - SA signals with NCBI Refseq genes
refseq.GRanges <- refseq %>%
  select(Symbol=name2, chrom, start=txStart, end=txEnd) %>%
  makeGRangesFromDataFrame(keep.extra.columns = TRUE)

SA.GRanges <- sentinels %>%
  mutate(Chrom=paste0("chr", chr)) %>%
  select(Chrom, start=pos, end=pos, sentinel) %>%
  makeGRangesFromDataFrame(keep.extra.columns = TRUE)

nearby_genes <- mergeByOverlaps(SA.GRanges, refseq.GRanges, maxgap=DISTANCE) %>%
  as_tibble %>%
  select(sentinel, Position=SA.GRanges.start, Symbol,
         txStart=refseq.GRanges.start, txEnd=refseq.GRanges.end,
         width=refseq.GRanges.width) %>%
  mutate(distance=ifelse(txStart <= Position & txEnd >= Position, 0,
                         ifelse(Position < txStart, txStart - Position, Position - txEnd))) %>%
  group_by(sentinel, Symbol) %>%
  filter(distance == min(distance)) %>%
  filter(width == max(width)) %>%
  slice(1) %>%
  ungroup

n_refseq <- refseq %>% getN(name2)
n_nearby_genes_sentinel <- nearby_genes %>% getN(sentinel)
n_nearby_genes_Symbol <- nearby_genes %>% getN(Symbol)

## Orthologs of mouse KO genes with asthma relevant phenotype - BROAD filter
#BROAD filter for top_level_mp_term == "respiratory system phenotype" give immunity and/or muscle related term as well:
ko_mouse_mp <- geno_pheno %>%
  filter(top_level_mp_term == "respiratory system phenotype" |
         top_level_mp_term == "immune system phenotype" |
         top_level_mp_term == "muscle phenotype") %>%
  select(top_level_mp_term, mp_term_name) %>%
  distinct %>% arrange(top_level_mp_term,mp_term_name)
#save the top level and suptypes level that I use:
fwrite(ko_mouse_mp,"/home/n/nnp5/PhD/PhD_project/Var_to_Gene/input/MP_TERM_approx_asthma_respimmmuscle.txt",quote=F, sep="\t")

#Actual filter for relevant top level and subtypes
ko_mouse <- geno_pheno %>%
  filter(top_level_mp_term == "respiratory system phenotype" |
         top_level_mp_term == "immune system phenotype" |
         top_level_mp_term == "muscle phenotype") %>%
  select(marker_symbol, mp_term_name) %>%
  distinct

ko_human <- orthologs %>%
  filter(support != "OrthoMCL") %>%
  select(mouse_symbol, human_symbol) %>%
  inner_join(ko_mouse, c("mouse_symbol" = "marker_symbol")) %>%
  mutate(overlap=human_symbol %in% nearby_genes$Symbol)

no_ortholog <- setdiff(ko_mouse$marker_symbol, ko_human$mouse_symbol)
ko_human_selected <- ko_human %>%
  filter(overlap)

## Hypergeometric test
hitInSample <- nrow(ko_human_selected)
hitInPop <- nrow(ko_human)
failInPop <- n_refseq - hitInPop
sampleSize <- n_nearby_genes_Symbol

p_resp_en <- phyper(q=hitInSample - 1,
                 m=hitInPop,
                 n=failInPop,
                 k=sampleSize,
                 lower.tail = FALSE)  #0.6866119

table3a.m <- matrix(c(hitInSample, hitInPop - hitInSample,
                      sampleSize - hitInSample, failInPop - sampleSize + hitInSample), 2, 2)
rownames(table3a.m) <- c("yes", "no")
colnames(table3a.m) <- c("asthma", "others")
table3a <- table3a.m %>%
  as_tibble(rownames = "Near severe asthma sentinel") %>%
  mutate(pc_asthma=(asthma/(asthma + others)) %>% percent)

pander(table3a, big.mark=",", caption=paste0("Table 3: For mouse knockout-causing genes, the proportion that are asthma related within ", DISTANCE/1000, "kb of an asthma function sentinel."))

### By-gene and by-SNP results
concat <- function(x) x[!is.na(x)] %>% unique %>% paste(collapse="; ")

implicated_mko <- right_join(nearby_genes, ko_human_selected, c("Symbol" = "human_symbol")) %>%
  add_column(MKO="MKO", .before = "txStart")

results <- left_join(sentinels, implicated_mko)

results_mko <- results %>%
  filter(!is.na(MKO)) %>%
  arrange(distance, chr, pos) %>%
  ungroup 

results_by_snp <- results %>%
  arrange(distance) %>%
  mutate_at(vars(MKO), ~ifelse(!is.na(.), .data$Symbol, NA)) %>%
  group_by(sentinel) %>%
  summarise_at(vars(MKO, distance), concat) %>%
  ungroup %>% 
  left_join(sentinels, .) 

results_by_snp_phyper <- results_mko %>%
  group_by(sentinel) %>%
  summarise(hitInSample_snp=n()) %>%
  ungroup %>%
  left_join(nearby_genes %>%
              group_by(sentinel) %>%
              summarise(sampleSize_snp=n())) %>%
  mutate(hyperG_p=phyper(q=hitInSample_snp - 1,
                 m=hitInPop,
                 n=failInPop,
                 k=sampleSize_snp,
                 lower.tail = FALSE)) %>%
  left_join(results_by_snp, .)

results_by_gene <- results %>%
  filter(!is.na(Symbol)) %>%
  arrange(distance) %>%
  mutate_at(vars(MKO), ~ifelse(!is.na(.), .data$sentinel, NA)) %>%
  group_by(Symbol) %>%
  summarise_at(vars(MKO, distance), concat) %>%
  ungroup 

results_by_snp_phyper %>% filter(!is.na(hyperG_p)) %>% count(hyperG_p < 0.05) %>%
  pander(caption="By SNP hypergeometric results")

out_base <- paste0("/scratch/gen1/nnp5/results_", DISTANCE/1000, "kb")
write_csv(results, paste0(out_base, ".csv"), na = "")
write_csv(results_mko, paste0(out_base, "_mko.csv"), na = "")
write_csv(results_by_snp_phyper, paste0(out_base, "_by_snp.csv"), na = "")
write_csv(results_by_gene, paste0(out_base, "_by_gene.csv"), na = "")
MP_TERM <- unique(ko_mouse_mp$mp_term_name)
inner_join(results_mko %>%
             select(-MKO, -txStart, -txEnd, -width, -overlap),
           results_by_snp_phyper %>% select(sentinel, hyperG_p)) %>%
  inner_join(all_geno_pheno %>%
               filter(mp_term_name %in% MP_TERM) %>%
               select(marker_symbol, mp_term_name, p_value, percentage_change, effect_size),
             c("mouse_symbol" = "marker_symbol")) %>%
  group_by(Symbol, mouse_symbol) %>%
  arrange(p_value) %>%
  slice(1) %>%
  write_csv("mouse_ko.csv")

#save genes only:
fwrite(as.data.frame(results_by_gene$Symbol), "/home/n/nnp5/PhD/PhD_project/Var_to_Gene/input/mouse_ko_genes_raw.txt", na = "",col.names=F,quote=F)