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
#MP_TERM <- "respiratory system phenotype"  

getN <- function(x, y) x %>% pull({{y}}) %>% n_distinct

orthologs <- read_tsv("human_mouse_hcop_fifteen_column.txt.gz")
all_geno_pheno <- read_csv("genotype-phenotype-assertions-ALL.csv.gz")
MP_TERM_lung <- unique(all_geno_pheno  %>% filter(grepl("lung",mp_term_name)) %>% select(mp_term_name))
MP_TERM_lung$MP_TERM <- "lung"
MP_TERM_airway <- unique(all_geno_pheno  %>% filter(grepl("airway",mp_term_name)) %>% select(mp_term_name))
MP_TERM_airway$MP_TERM <- "airway"
MP_TERM_muscle <- unique(all_geno_pheno  %>% filter(grepl("muscle",mp_term_name)) %>% select(mp_term_name))
MP_TERM_muscle$MP_TERM <- "muscle"
MP_TERM_imm <- unique(all_geno_pheno  %>% filter(grepl("imm",mp_term_name)) %>% select(mp_term_name))
MP_TERM_imm$MP_TERM <- "imm"
MP_TERM_epith <-  unique(all_geno_pheno  %>% filter(grepl("epith",mp_term_name)) %>% select(mp_term_name))
MP_TERM_epith$MP_TERM <- "epith"
MP_TERM_bronchoconstri <-  unique(all_geno_pheno  %>% filter(grepl("bronchoconstri",mp_term_name)) %>% select(mp_term_name))
MP_TERM_bronchoconstri$MP_TERM <- "bronchoconstri"
MP_TERM_pulm <-  unique(all_geno_pheno  %>% filter(grepl("pulm",mp_term_name)) %>% select(mp_term_name))
MP_TERM_pulm$MP_TERM <- "pulm"

#put all data frames into list
df_list <- list(MP_TERM_airway,MP_TERM_epith,MP_TERM_muscle,MP_TERM_bronchoconstri,MP_TERM_imm,MP_TERM_lung,MP_TERM_pulm)
#merge all data frames in list
MP_TERM <- df_list %>% reduce(full_join, by = join_by(mp_term_name, MP_TERM))
fwrite(MP_TERM,"/home/n/nnp5/PhD/PhD_project/Var_to_Gene/input/MP_TERM_approx_asthma.txt",quote=F, sep="\t")

sentinels <- read_tsv("/scratch/gen1/jc824/TSH/novel_signals/TSH_signal_list.txt")
setnames(sentinels,"MarkerName","sentinel")
setnames(sentinels,"Chromosome","chr")
setnames(sentinels,"Position","pos")

ucsc <- browserSession("UCSC")
genome(ucsc) <- "hg19"
refseq <- ucsc %>%
  ucscTableQuery(track="NCBI RefSeq", table="refGene") %>%
  getTable %>%
  as_tibble

## Mouse genotype-phenotype data
geno_pheno <- all_geno_pheno %>%
  separate(top_level_mp_term_name, c("top_level_mp_term", "top_level_mp_subtype"), ",")

geno_pheno %>% 
  group_by(top_level_mp_term) %>%
  summarise(nlines=n(), ngenes=n_distinct(marker_symbol)) %>%
  drop_na %>%
  pander(big.mark=",")

## Overlap of TSH signals with NCBI Refseq genes
refseq.GRanges <- refseq %>%
  select(Symbol=name2, chrom, start=txStart, end=txEnd) %>%
  makeGRangesFromDataFrame(keep.extra.columns = TRUE)

TSH.GRanges <- sentinels %>%
  mutate(Chrom=paste0("chr", chr)) %>%
  select(Chrom, start=pos, end=pos, sentinel) %>%
  makeGRangesFromDataFrame(keep.extra.columns = TRUE)

nearby_genes <- mergeByOverlaps(TSH.GRanges, refseq.GRanges, maxgap=DISTANCE) %>%
  as_tibble %>%
  select(sentinel, Position=TSH.GRanges.start, Symbol,
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

## Orthologs of mouse KO genes with thyro relevant phenotype
ko_mouse <- all_geno_pheno %>%
  filter(mp_term_name %in% MP_TERM) %>%
  select(marker_symbol,mp_term_name) %>%
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
colnames(table3a.m) <- c("thyro", "others")
table3a <- table3a.m %>%
  as_tibble(rownames = "Near TSH sentinel") %>%
  mutate(pc_thyro=(thyro/(thyro + others)) %>% percent)

pander(table3a, big.mark=",", caption=paste("Table 3: For mouse knockout-causing genes, the proportion that are thyro related within $\\pm$", DISTANCE/1000, "kb of a lung function sentinel."))

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

out_base <- paste0("results_", DISTANCE/1000, "kb")
write_csv(results, paste0(out_base, ".csv"), na = "")
write_csv(results_mko, paste0(out_base, "_mko.csv"), na = "")
write_csv(results_by_snp_phyper, paste0(out_base, "_by_snp.csv"), na = "")
write_csv(results_by_gene, paste0(out_base, "_by_gene.csv"), na = "")
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
  