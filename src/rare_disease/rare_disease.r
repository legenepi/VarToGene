library(GenomicRanges)
library(rtracklayer)
library(tidyverse)
library(pander)
library(scales)
library(readxl)
library(data.table)
# R/4.0.0

setwd("/data/gen1/TSH/rare_disease/")
# Define distance to look up nearby rare disease associated genes
DISTANCE <- 5e5
# funtion to work out the number unique element of the yth column of data table x
getN <- function(x, y) x %>% pull({{y}}) %>% n_distinct

## retrieve genomic reference for genes
ucsc = browserSession("UCSC")  # methods for getting browser sessions
genome(ucsc) <- "hg19"  # setting the sequence information stored in ucsc
# genomic reference for genes
refseq <- ucsc %>%
  ucscTableQuery(track="NCBI RefSeq", table="refGene") %>%
  getTable %>%
  as_tibble

## genes and rare disease association from orphadata 
# http://www.orphadata.org/data/xml/en_product6.xml
orphanet_genes <- read_xlsx("en_product6.xlsx") %>%
  select(Symbol, Gene=Name12, OrphaCode, Disease=Name, SourceOfValidation,
         DisorderGeneAssociationType=Name24) %>%
  distinct
hpo_f_levels <- c("Obligate (100%)","Very frequent (99-80%)","Frequent (79-30%)","Occasional (29-5%)","Very rare (<4-1%)","Excluded (0%)")
# hpo: human phenotype ontology provides a standardized vocabulary of phenotypic abnormalities encounterd in human disease
# http://www.orphadata.org/data/xml/en_product4.xml
orphanet_hpo <- read_xlsx("en_product4.xlsx") %>%
  select(OrphaCode, HPOId, HPOTerm, HPOFrequency=Name15) %>%
  distinct %>%
  mutate(HPOFrequency=factor(HPOFrequency, levels=hpo_f_levels))

## TSH signals
sentinels <- read_tsv("/data/gen1/TSH/novel_signals/TSH_signal_list.txt")
setnames(sentinels,"MarkerName","sentinel")
setnames(sentinels,"Chromosome","chr")
setnames(sentinels,"Position","pos")

## Overlap of TSH function signals with NCBI Refseq genes
# Create a GRanges object for reference of genes
refseq.GRanges <- refseq %>%
  mutate(chrom=sub("^chr", "", chrom)) %>%
  select(Symbol=name2, chrom, start=txStart, end=txEnd) %>%
  makeGRangesFromDataFrame(keep.extra.columns = TRUE) # takes a data-frame-like object as input and tries to automatically find the columns that describe genomic ranges. It returns them as a GRanges object 
# Create a GRanges object for our signals
sentinels.GRanges <- sentinels %>%
  select(chr, start=pos, sentinel) %>%
  makeGRangesFromDataFrame(end.field = "start", keep.extra.columns = TRUE)
# overlap checking
nearby_genes <- mergeByOverlaps(sentinels.GRanges, refseq.GRanges, maxgap=DISTANCE) %>%
  as_tibble %>%
  select(sentinel, Position=sentinels.GRanges.start, Symbol,
         txStart=refseq.GRanges.start, txEnd=refseq.GRanges.end,
         width=refseq.GRanges.width) %>%
  mutate(distance=ifelse(txStart <= Position & txEnd >= Position, 0,
                         ifelse(Position < txStart, txStart - Position, Position - txEnd))) %>%
  group_by(sentinel, Symbol) %>%   # most data operation are done one groups defined by variables
  filter(distance == min(distance)) %>%
  filter(width == max(width)) %>%
  slice(1) %>%
  ungroup

n_refseq <- refseq %>% getN(name2)
n_nearby_genes_Sentinel <- nearby_genes %>% getN(sentinel)
n_nearby_genes_Symbol <- nearby_genes %>% getN(Symbol)
n_genes <- orphanet_genes %>% getN(Symbol)
n_disease <- orphanet_genes %>% getN(OrphaCode)
n_disease_with_hpo <- orphanet_hpo %>% getN(OrphaCode)
n_hpo <- orphanet_hpo %>% getN(HPOId)

# Merge Orphanet genes with HPO terms and keep only Disease-causing germline genes
orphanet <- left_join(orphanet_genes, orphanet_hpo) %>%
  filter(grepl("Disease-causing germline", DisorderGeneAssociationType)) %>%
  mutate(HPOFrequency = fct_explicit_na(HPOFrequency, "No HPO term"))

n_gene_disease_c <- orphanet %>% group_by(Symbol, Disease) %>% n_groups 
n_genes_c <- orphanet %>% getN(Symbol)
n_disease_c <- orphanet %>% getN(Disease)
n_disease_hpo_c <- orphanet %>% group_by(Disease, HPOId) %>% n_groups 
n_disease_with_hpo_c <- orphanet %>% filter(!is.na(HPOId)) %>% getN(Disease)
n_hpo_c <- orphanet %>% getN(HPOId)

orphanet %>%
  count(HPOFrequency) %>%
  mutate(pc=(n/sum(n)) %>% percent) %>%
  pander(big.mark=",", caption=paste("Table 1a: Distribution of HPO term frequencies across",
  nrow(orphanet) %>% comma, "rows of gene-disease-hpo associations for Disease-causing germline genes"))

orphanet %>%
  group_by(Symbol) %>%
  arrange(HPOFrequency) %>%
  slice(1) %>%
  ungroup %>%
  count(HPOFrequency) %>%
  mutate(pc=(n/sum(n)) %>% percent) %>%
  pander(big.mark=",", caption=paste("Table 1b: Distribution of highest HPO term frequency per gene across",
                                     n_genes_c %>% comma, "Disease-causing germline genes"))

## Filtering for thyroid terms
red <- hpo_f_levels[1:2]
yellow <- hpo_f_levels[3]
green <- hpo_f_levels[4:5]

orphanet <- orphanet %>%
  mutate(thyro_HPO=grepl("thyro", HPOTerm, ignore.case = TRUE) &
                   HPOFrequency != "Excluded (0%)",
         thyro_disease=grepl("thyro", Disease, ignore.case = TRUE),
         parathyro=!grepl("parathyro", HPOTerm, ignore.case = TRUE) &
                   !grepl("parathyro", Disease, ignore.case = TRUE),
         thyro=(thyro_disease | thyro_HPO) & parathyro,
         overlap=Symbol %in% nearby_genes$Symbol,
         evidence=ifelse(!thyro, "Not thyroid related",
                         ifelse(thyro_disease | HPOFrequency %in% red,
                         "Disease name, obligate or very frequent",
                         ifelse(HPOFrequency %in% yellow, "Frequent",
                                ifelse(HPOFrequency %in% green, "Occasional or rare", "Excluded")))) %>%
           factor(levels = c("Disease name, obligate or very frequent",
                             "Frequent",
                             "Occasional or rare",
                             "Excluded",
                             "Not thyroid related"))) 

n_thyro_genes_disease <- orphanet %>% filter(thyro_disease) %>% getN(Symbol)
n_thyro_genes_hpo <- orphanet %>% filter(thyro_HPO) %>% getN(Symbol)

orphanet_thyro <- orphanet %>%
  filter(thyro)

n_thyro_genes <- orphanet_thyro %>% getN(Symbol)

orphanet_thyro_selected <- orphanet_thyro %>%
  filter(overlap) %>%
  select(-overlap)

n_thyro_genes_selected <- orphanet_thyro_selected %>% getN(Symbol)

table2a <- orphanet %>%
  group_by(Symbol) %>%
  summarise_at(vars(thyro, overlap), any) %>%
  ungroup %>%
  count(thyro, overlap) %>%
  pivot_wider(names_from = "overlap", values_from = "n", names_prefix = "overlap_") %>%
  mutate_at(vars(starts_with("overlap")), list(pc=~(./sum(.)) %>% percent))

table2a_chisq <- table2a %>% select(2:3) %>% chisq.test

table2b <- orphanet %>%
  group_by(Symbol) %>%
  arrange(evidence) %>%
  slice(1) %>%
  ungroup %>%
  count(evidence, overlap) %>%
  pivot_wider(names_from = "overlap", values_from = "n", names_prefix = "overlap_",
              values_fill = list(n=0)) %>%
  mutate_at(vars(starts_with("overlap")), list(pc=~(./sum(.)) %>% percent))

table2b_chisq <- table2b %>% select(2:3) %>% chisq.test

pander(table2a, big.mark=",",
       caption=paste("Table 2a: For", n_genes_c %>% comma, "rare Disease-causing, germline genes, the numbers of thyroid and non-thyroid related genes vs. numbers overlapping TSH signals $\\pm$", DISTANCE/1000, "kb. $\\chi^2$ P =", table2a_chisq$p.value %>% round(3)))

pander(table2b, big.mark=",",
       caption=paste("Table 2b: For", n_genes_c %>% comma, "rare Disease-causing, germline genes, the numbers of genes in each thyroid evidence category vs. numbers overlapping TSH signals $\\pm$", DISTANCE/1000, "kb. $\\chi^2$ P =", table2b_chisq$p.value %>% round(3)), split.table=Inf)

## Hypergeometric test
hitInSample <- n_thyro_genes_selected
hitInPop <- n_thyro_genes
failInPop <- n_refseq - n_thyro_genes
sampleSize <- n_nearby_genes_Symbol

p_thyro_en <- phyper(q=hitInSample - 1,
                 m=hitInPop,
                 n=failInPop,
                 k=sampleSize,
                 lower.tail = FALSE)  # P value for hypergeometric test: 0.001188

table3a.m <- matrix(c(hitInSample, hitInPop - hitInSample,
                      sampleSize - hitInSample, failInPop - sampleSize + hitInSample), 2, 2)
rownames(table3a.m) <- c("yes", "no")
colnames(table3a.m) <- c("thyroid", "others")
table3a <- table3a.m %>%
  as_tibble(rownames = "Near TSH sentinel") %>%
  mutate(pc_thyro=(thyroid/(thyroid + others)) %>% percent)
pander(table3a, big.mark=",", caption=paste("Table 3: For rare disease-causing genes, the proportion that are thyroid related within $\\pm$", DISTANCE/1000, "kb of a TSH sentinel."))

## By-gene and by-SNP results
concat <- function(x) {
  x <- x[!is.na(x)]
  if (length(x) > 0 & is.numeric(x)) {
    min(x, na.rm=T) %>% paste
  } else {
    x %>% unique %>% paste(collapse="; ")
  }
}

orphanet_thyro_selected_wide <- orphanet_thyro_selected %>%
  mutate_at(vars(starts_with("HPO")), ~ifelse(.data$thyro_HPO, ., NA)) %>%
  arrange(evidence) %>%
  group_by(Symbol, Gene) %>%
  summarise(Evidence=evidence[1],
            nDisease=n_distinct(Disease, na.rm = TRUE),
            Diseases=concat(Disease),
            Validation=concat(SourceOfValidation),
            nHPOTerm=n_distinct(HPOTerm, na.rm = TRUE),
            HPOTerms=concat(HPOTerm)) %>%
  ungroup %>%
  arrange(Evidence)
  implicated_rare <- right_join(nearby_genes, orphanet_thyro_selected_wide) %>%
  add_column(Rare="Rare", .before = "txStart")

results <- implicated_rare %>% select(-Position) %>% left_join(sentinels, .)

results_rare <- results %>% filter(!is.na(Rare)) %>% arrange(Evidence, distance, chr, pos) %>% ungroup 

results_by_snp <- results %>%
  arrange(Evidence, distance) %>%
  mutate_at(vars(Rare), ~ifelse(!is.na(.), .data$Symbol, NA)) %>%
  group_by(sentinel) %>%
  summarise_at(vars(Rare, Evidence, distance), concat) %>%
  ungroup %>% 
  left_join(sentinels, .)

results_by_snp_phyper <- results_rare %>%
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
  arrange(Evidence, distance) %>%
  mutate_at(vars(Rare), ~ifelse(!is.na(.), .data$sentinel, NA)) %>%
  group_by(Symbol) %>%
  summarise_at(vars(Rare, Evidence, distance, Diseases, HPOTerms), concat)

out_base <- paste0("results_", DISTANCE/1000, "kb")
write_csv(results, paste0(out_base, ".csv"), na = "")
write_csv(results_rare, paste0(out_base, "_rare.csv"), na = "")
write_csv(results_by_snp_phyper, paste0(out_base, "_by_snp.csv"), na = "")
write_csv(results_by_gene, paste0(out_base, "_by_gene.csv"), na = "")