#!/usr/bin/env Rscript

#Rationale: Find nearby rare mendelian genes that are within +/- 500Kb from SA GWAS sentinels.

library(GenomicRanges)
library(rtracklayer)
library(tidyverse)
library(pander)
library(scales)
library(readxl)
library(data.table)

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
#https://www.orphadata.com/data/xml/en_product6.xml
orphanet_genes <- read_xlsx("/home/n/nnp5/PhD/PhD_project/Var_to_Gene/input/en_product6.xlsx") %>%
  select(Symbol, Gene=Name12, OrphaCode, Disease=Name, SourceOfValidation,
         DisorderGeneAssociationType=Name24) %>%
  distinct

# hpo: human phenotype ontology provides a standardized vocabulary of phenotypic abnormalities encountered in human disease
#https://www.orphadata.org/data/xml/en_product4.xml
#create the levels for HPOFrequency column:
hpo_f_levels <- c("Obligate (100%)","Very frequent (99-80%)","Frequent (79-30%)","Occasional (29-5%)","Very rare (<4-1%)","Excluded (0%)")
orphanet_hpo <- read_xlsx("/home/n/nnp5/PhD/PhD_project/Var_to_Gene/input/en_product4.xlsx") %>%
  select(OrphaCode, HPOId, HPOTerm, HPOFrequency=Name15) %>%
  distinct %>%
  mutate(HPOFrequency=factor(HPOFrequency, levels=hpo_f_levels))

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

## Overlap of GWAS function signals with NCBI Refseq genes
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

## Filtering for asthma terms
red <- hpo_f_levels[1:2]
yellow <- hpo_f_levels[3]
green <- hpo_f_levels[4:5]

#HPOId term for asthma: HP:0002099
orphanet <- orphanet %>%
  mutate(asthma_HPO=grepl("asthma", HPOTerm, ignore.case = TRUE) &
                   HPOFrequency != "Excluded (0%)",
         asthma=asthma_HPO,
         overlap=Symbol %in% nearby_genes$Symbol,
         evidence=ifelse(!asthma, "Not asthma related",
                         ifelse(HPOFrequency %in% red,
                         "Disease name, obligate or very frequent",
                         ifelse(HPOFrequency %in% yellow, "Frequent",
                                ifelse(HPOFrequency %in% green, "Occasional or rare", "Excluded")))) %>%
           factor(levels = c("Disease name, obligate or very frequent",
                             "Frequent",
                             "Occasional or rare",
                             "Excluded",
                             "Not asthma related")))

#n_asthma_genes_disease <- orphanet %>% filter(asthma_disease) %>% getN(Symbol) #not available for asthma
n_asthma_genes_hpo <- orphanet %>% filter(asthma_HPO) %>% getN(Symbol)

orphanet_asthma <- orphanet %>%
  filter(asthma)

n_asthma_genes <- orphanet_asthma %>% getN(Symbol)

orphanet_asthma_selected <- orphanet_asthma %>%
  filter(overlap) %>%
  select(-overlap)

n_asthma_genes_selected <- orphanet_asthma_selected %>% getN(Symbol)

table2a <- orphanet %>%
  group_by(Symbol) %>%
  summarise_at(vars(asthma, overlap), any) %>%
  ungroup %>%
  count(asthma, overlap) %>%
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
       caption=paste0("Table 2a: For", n_genes_c %>% comma, "rare Disease-causing, germline genes, the numbers of asthma and non-asthma related genes vs. numbers overlapping severe asthma signals ", DISTANCE/1000, "kb. chi^2 P =", table2a_chisq$p.value %>% round(3)))

pander(table2b, big.mark=",",
       caption=paste0("Table 2b: For", n_genes_c %>% comma, "rare Disease-causing, germline genes, the numbers of genes in each asthma evidence category vs. numbers overlapping severe asthma signals ", DISTANCE/1000, "kb. chi^2 P =", table2b_chisq$p.value %>% round(3)), split.table=Inf)

## Hypergeometric test
hitInSample <- n_asthma_genes_selected
hitInPop <- n_asthma_genes
failInPop <- n_refseq - n_asthma_genes
sampleSize <- n_nearby_genes_Symbol

p_asthma_en <- phyper(q=hitInSample - 1,
                 m=hitInPop,
                 n=failInPop,
                 k=sampleSize,
                 lower.tail = FALSE)  # P value for hypergeometric test: 0.001188

table3a.m <- matrix(c(hitInSample, hitInPop - hitInSample,
                      sampleSize - hitInSample, failInPop - sampleSize + hitInSample), 2, 2)
rownames(table3a.m) <- c("yes", "no")
colnames(table3a.m) <- c("asthma", "others")
table3a <- table3a.m %>%
  as_tibble(rownames = "Near severe asthma sentinel") %>%
  mutate(pc_asthma=(asthma/(asthma + others)) %>% percent)
pander(table3a, big.mark=",", caption=paste0("Table 3: For rare disease-causing genes, the proportion that are asthma related within ", DISTANCE/1000, "kb of a severe asthma sentinel."))

## By-gene and by-SNP results
concat <- function(x) {
  x <- x[!is.na(x)]
  if (length(x) > 0 & is.numeric(x)) {
    min(x, na.rm=T) %>% paste
  } else {
    x %>% unique %>% paste(collapse="; ")
  }
}

orphanet_asthma_selected_wide <- orphanet_asthma_selected %>%
  mutate_at(vars(starts_with("HPO")), ~ifelse(.data$asthma_HPO, ., NA)) %>%
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
  implicated_rare <- right_join(nearby_genes, orphanet_asthma_selected_wide) %>%
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

out_base <- paste0("/scratch/gen1/nnp5/Var_to_Gen_tmp/rare_disease/results_", DISTANCE/1000, "kb")
write_csv(results, paste0(out_base, ".csv"), na = "")
write_csv(results_rare, paste0(out_base, "_rare.csv"), na = "")
write_csv(results_by_snp_phyper, paste0(out_base, "_by_snp.csv"), na = "")
write_csv(results_by_gene, paste0(out_base, "_by_gene.csv"), na = "")
#save genes only:
fwrite(as.data.frame(results_by_gene$Symbol), "/home/n/nnp5/PhD/PhD_project/Var_to_Gene/input/rare_disease_genes_raw.txt", na = "",col.names=F,quote=F)