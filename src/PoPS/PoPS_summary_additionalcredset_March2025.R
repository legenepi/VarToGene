#!/usr/bin/env Rscript

#Rationale: Look at PoPS results and create a summary of them - additional credset vars March 2025
##We currently suggest taking the highest scoring gene in each GWAS locus.
##You could further filter this set to only include genes in top 10% of PoP scores across all genes.
##Negative scores generally mean low evidence -- this is the predicted MAGMA z-score!!’ (https://github.com/FinucaneLab/pops/issues/4)
##From TSH paper: ‘prioritized genes for all autosomal TSH sentinel variants within a 500kb (±250kb) window of
#the sentinel variant and reported the top prioritised genes in the region.
#If there was no gene prioritized within a 500kb window of the sentinel,
#we reported any top prioritized genes within a 1Mb window (Supplementary Data 19).’.


library(data.table)
library(tidyverse)

#Create file with PIP sentinel variants from the Credible_sets.xlsx file and import it:
sig_list <- fread("input/Additional_credset_snps_March2025/Locus_PIPSentinel_chr_pos.txt") %>% rename(sentinel = PIP_Sentinel, locus = 'Replicated locus')
sig_list$chr <- as.numeric(sig_list$chr)
sig_list$pos <- as.numeric(sig_list$pos)
sig_list <- sig_list %>% mutate(locus = paste0(sig_list$locus,"_", sig_list$sentinel))
locus <- sig_list$locus
gene_loc <- fread("/data/gen1/LF_HRC_transethnic/PoPS/data/gene_loc.txt")
gene_features <- fread("/data/gen1/LF_HRC_transethnic/PoPS/data/PoPS.features.txt.gz")


##### sentinel 500 kb total window #####
result <- data.frame(matrix(ncol = 16,nrow = 0))
for(i in locus){
    locus_sig_list <- sig_list %>% filter(locus == as.character(i))
    sentinel <- locus_sig_list %>% filter(locus == as.character(i)) %>% select(sentinel)
    chr <- locus_sig_list %>% filter(locus == as.character(i)) %>% select(chr)
    pos <- locus_sig_list %>% filter(locus == as.character(i)) %>% select(pos)
    print(sentinel)
    print(chr)
    print(pos)
    trait <- "SA"
    if(chr!="X"){
        chr <- as.integer(chr)
        pos <- as.integer(pos)
        pos1 <- pos-250000 # window set
        pos2 <- pos+250000 # window set
        if(pos1<0){
            pos1 <- 0
        }
        genes <- gene_loc[gene_loc$CHR==chr,]
        genes <- genes[(genes$START>=pos1&genes$START<=pos2)|(genes$END>=pos1&genes$END<=pos2),]
        n_genes <- nrow(genes)
        if(n_genes!=0){ 
            result1 <- data.frame(matrix(ncol = 16,nrow = 0))          
            for(j in 1:n_genes){
                result2 <- data.frame(matrix(ncol = 16,nrow = 0))
                gene <- genes$ENSGID[j]
                result2[1,1] <- gene
                result2[1,2] <- sentinel
                PoPS <- fread(paste0("/scratch/gen1/nnp5/Var_to_Gen_tmp/pops/",trait,".",chr,".results"))
                score <- PoPS[PoPS$ENSGID==gene,]$Score[1]
                result2[1,3] <- score
                result2[1,4] <- i
                beta <- fread(paste0("/scratch/gen1/nnp5/Var_to_Gen_tmp/pops/",trait,".",chr,".coefs"))
                X_all <- gene_features[ENSGID==gene,]
                X <- select(X_all,beta$Feature)
                X <- as.data.table(t(X))
                beta$X <- X$V1
                beta <- beta[,contribute:=X*beta_hat]
                beta <- beta[order(abs(contribute),decreasing=TRUE)]
                result2[1,7] <- beta$Feature[1]
                result2[1,8] <- beta$Feature[2]
                result2[1,9] <- beta$Feature[3]
                result2[1,10] <- beta$Feature[4]
                result2[1,11] <- beta$Feature[5]
                result2[1,12] <- beta$Feature[6]
                result2[1,13] <- beta$Feature[7]
                result2[1,14] <- beta$Feature[8]
                result2[1,15] <- beta$Feature[9]
                result2[1,16] <- beta$Feature[10]
                result1 <<- rbind(result1,result2)      
            }
            result1 <- as.data.table(result1)
            result1 <- result1[order(X3,decreasing=TRUE)]
            result1$X6 <- FALSE
            result1$X6[1] <- TRUE
            for(j in 1:n_genes){
                result1$X5[j] <- j
            }
            result <<- rbind(result,result1)
        }
    }
}
names(result) <- c("ENSGID","sentinel","PoPS_score","signal_id","gene_rank","prioritized","Feature1","Feature2","Feature3","Feature4","Feature5","Feature6","Feature7","Feature8","Feature9","Feature10")
result <- result[order(signal_id,gene_rank)]
write.table(result,"/scratch/gen1/nnp5/Var_to_Gen_tmp/pops/results/additional_credsetMarch2025_all_result_summary_250kb_window.txt", row.names=F, quote=F, sep="\t")
#save results for only prioritised genes:
write.table(result[gene_rank==1,],"/scratch/gen1/nnp5/Var_to_Gen_tmp/pops/results/additional_credsetMarch2025_prioritized_genes_250kb_window.txt", row.names=F, quote=F, sep="\t")
#save prioritised genes only:
genes_250 <- unique(result %>% filter(gene_rank == 1) %>% select(ENSGID))
write.table(genes_250,"/scratch/gen1/nnp5/Var_to_Gen_tmp/pops/results/additional_credsetMarch2025_pops_var2genes_raw_250kbwindow.txt", row.names=F, quote=F, sep="\t", col.names=F)

##### sentinel 1 mb total window #####
result <- data.frame(matrix(ncol = 16,nrow = 0))
for(i in locus){
    locus_sig_list <- sig_list %>% filter(locus == as.character(i))
    sentinel <- locus_sig_list %>% filter(locus == as.character(i)) %>% select(sentinel)
    chr <- locus_sig_list %>% filter(locus == as.character(i)) %>% select(chr)
    pos <- locus_sig_list %>% filter(locus == as.character(i)) %>% select(pos)
    print(sentinel)
    print(chr)
    print(pos)
    trait <- "SA"
    if(chr!="X"){
        chr <- as.integer(chr)
        pos <- as.integer(pos)
        pos1 <- pos-500000 # window set
        pos2 <- pos+500000 # window set
        if(pos1<0){
            pos1 <- 0
        }
        genes <- gene_loc[gene_loc$CHR==chr,]
        genes <- genes[(genes$START>=pos1&genes$START<=pos2)|(genes$END>=pos1&genes$END<=pos2),]
        n_genes <- nrow(genes)
        if(n_genes!=0){
            result1 <- data.frame(matrix(ncol = 16,nrow = 0))
            for(j in 1:n_genes){
                result2 <- data.frame(matrix(ncol = 16,nrow = 0))
                gene <- genes$ENSGID[j]
                result2[1,1] <- gene
                result2[1,2] <- sentinel
                PoPS <- fread(paste0("/scratch/gen1/nnp5/Var_to_Gen_tmp/pops/",trait,".",chr,".results"))
                score <- PoPS[PoPS$ENSGID==gene,]$Score[1]
                result2[1,3] <- score
                result2[1,4] <- i
                beta <- fread(paste0("/scratch/gen1/nnp5/Var_to_Gen_tmp/pops/",trait,".",chr,".coefs"))
                X_all <- gene_features[ENSGID==gene,]
                X <- select(X_all,beta$Feature)
                X <- as.data.table(t(X))
                beta$X <- X$V1
                beta <- beta[,contribute:=X*beta_hat]
                beta <- beta[order(abs(contribute),decreasing=TRUE)]
                result2[1,7] <- beta$Feature[1]
                result2[1,8] <- beta$Feature[2]
                result2[1,9] <- beta$Feature[3]
                result2[1,10] <- beta$Feature[4]
                result2[1,11] <- beta$Feature[5]
                result2[1,12] <- beta$Feature[6]
                result2[1,13] <- beta$Feature[7]
                result2[1,14] <- beta$Feature[8]
                result2[1,15] <- beta$Feature[9]
                result2[1,16] <- beta$Feature[10]
                result1 <<- rbind(result1,result2)
            }
            result1 <- as.data.table(result1)
            result1 <- result1[order(X3,decreasing=TRUE)]
            result1$X6 <- FALSE
            result1$X6[1] <- TRUE
            for(j in 1:n_genes){
                result1$X5[j] <- j
            }
            result <<- rbind(result,result1)
        }
    }
}
names(result) <- c("ENSGID","sentinel","PoPS_score","signal_id","gene_rank","prioritized","Feature1","Feature2","Feature3","Feature4","Feature5","Feature6","Feature7","Feature8","Feature9","Feature10")
result <- result[order(signal_id,gene_rank)]
write.table(result,"/scratch/gen1/nnp5/Var_to_Gen_tmp/pops/results/additional_credsetMarch2025_all_result_summary_500kb_window.txt", row.names=F, quote=F, sep="\t")
write.table(result[gene_rank==1,],"/scratch/gen1/nnp5/Var_to_Gen_tmp/pops/results/additional_credsetMarch2025_prioritized_genes_500kb_window.txt", row.names=F, quote=F, sep="\t")
#save prioritised genes only:
genes_500 <- unique(result %>% filter(gene_rank == 1) %>% select(ENSGID))
write.table(genes_500,"/scratch/gen1/nnp5/Var_to_Gen_tmp/pops/results/additional_credsetMarch2025_pops_var2genes_raw_500kbwindow.txt", row.names=F, quote=F, sep="\t", col.names=F)

################# mapping features ###############
# Get gene information using biomaRt for all genes
library(biomaRt)         # Requires R/4.1.0
#250Kb window results
pops <- fread("/scratch/gen1/nnp5/Var_to_Gen_tmp/pops/results/additional_credsetMarch2025_prioritized_genes_250kb_window.txt")
ensembl <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl", version = "GRCh37")
genes <- as_tibble(getBM(attributes=c('hgnc_symbol', 'ensembl_gene_id', 'strand'),
  filters=c('ensembl_gene_id'),
  values=list(sort(unique(pops$ENSGID))),
  mart=ensembl)) %>%
  rename(gene_symbol = hgnc_symbol, gene_strand = strand)

genes <- as.data.table(genes)
setnames(genes,"ensembl_gene_id","ENSGID")
merged <- merge(pops,genes,by="ENSGID")
write.table(merged,"/scratch/gen1/nnp5/Var_to_Gen_tmp/pops/results/additional_credsetMarch2025_all_result_summary_250kb_window_gene_mapped.txt", row.names=F, quote=F, sep="\t")

merged <- as.data.frame(merged[PoPS_score>0,])
write.table(merged,"/scratch/gen1/nnp5/Var_to_Gen_tmp/pops/results/additional_credsetMarch2025_all_results_merged_table.txt", row.names=F, quote=F, sep="\t")

#save final list of genes:
genes <- unique(merged$gene_symbol)
write.table(genes,"/home/n/nnp5/PhD/PhD_project/Var_to_Gene/input/Additional_credset_snps_March2025/pops_var2genes_raw_additional_credsetMarch2025.txt", row.names=F, quote=F, sep="\t")
