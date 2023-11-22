#!/usr/bin/env Rscript

#Rationale: Look at PoPS results and create a summary of them
##‘We currently suggest taking the highest scoring gene in each GWAS locus.
##You could further filter this set to only include genes in top 10% of PoP scores across all genes.
##Negative scores generally mean low evidence -- this is the predicted MAGMA z-score!!’ (https://github.com/FinucaneLab/pops/issues/4)
##From TSH paper: ‘prioritized genes for all autosomal TSH sentinel variants within a 500kb (±250kb) window of
#the sentinel variant and reported the top prioritised genes in the region.
#If there was no gene prioritized within a 500kb window of the sentinel,
#we reported any top prioritized genes within a 1Mb window (Supplementary Data 19).’.

library(data.table)
library(tidyverse)

sig_list <- fread("/data/gen1/UKBiobank_500K/severe_asthma/Noemi_PhD/data/replsugg_valid_credset.txt")
setnames(sig_list,"chromosome","chr")
setnames(sig_list,"position","pos")
sig_list$sentinel <- paste0(sig_list$chr,"_",sig_list$pos,"_",sig_list$allele1,"_",sig_list$allele2)
locus <- unique(sig_list$locus)
gene_loc <- fread("/data/gen1/LF_HRC_transethnic/PoPS/data/gene_loc.txt")
gene_features <- fread("/data/gen1/LF_HRC_transethnic/PoPS/data/PoPS.features.txt.gz")

##### sentinel 500 kb total window #####
result <- data.frame(matrix(ncol = 16,nrow = 0))
for(i in locus){
    locus_sig_list <- sig_list %>% filter(locus == as.character(i))
    sentinel <- locus_sig_list %>% filter(PIP_average == max(locus_sig_list$PIP_average)) %>% select(sentinel)
    chr <- locus_sig_list %>% filter(PIP_average == max(locus_sig_list$PIP_average)) %>% select(chr)
    pos <- locus_sig_list %>% filter(PIP_average == max(locus_sig_list$PIP_average)) %>% select(pos)
    print(sentinel)
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
                PoPS <- fread(paste0("/scratch/gen1/nnp5/Var_to_Gen_tmp/",trait,".",chr,".results"))
                score <- PoPS[PoPS$ENSGID==gene,]$Score[1]
                result2[1,3] <- score
                result2[1,4] <- i
                beta <- fread(paste0("/scratch/gen1/nnp5/Var_to_Gen_tmp/",trait,".",chr,".coefs"))
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
LOOK AT THIS CHUNK: MAYBE NOT NEEDED TO BE DONE:
#map1 <- fread("/scratch/gen1/nrgs1/ensGene.txt.gz")
#map2 <- fread("/scratch/gen1/nrgs1/ensemblToGeneName.txt.gz")
#map <- merge(map1,map2,by="name")
#map <- select(map,"name2","value")
#names(map) <- c("ENSGID","value")
#map <- unique(map)
#result <- unique(result)
#result_all <- merge(result,map,by="ENSGID",all.x=T,all.y=F)
#result_all <- result_all[order(signal_id,gene_rank)]
#write.table(result_all,"/scratch/gen1/nnp5/Var_to_Gen_tmp/pops/results/result_summary_250kb_window.txt", row.names=F, quote=F, sep="\t")
write.table(result,"/scratch/gen1/nnp5/Var_to_Gen_tmp/pops/results/all_result_summary_250kb_window.txt", row.names=F, quote=F, sep="\t")
write.table(result[gene_rank==1,],"/scratch/gen1/nnp5/Var_to_Gen_tmp/pops/results/prioritized_genes_250kb_window.txt", row.names=F, quote=F, sep="\t")

##### sentinel 1 mb total window #####
result <- data.frame(matrix(ncol = 16,nrow = 0))
for(i in locus){
    locus_sig_list <- sig_list %>% filter(locus == as.character(i))
    sentinel <- locus_sig_list %>% filter(PIP_average == max(locus_sig_list$PIP_average)) %>% select(sentinel)
    chr <- locus_sig_list %>% filter(PIP_average == max(locus_sig_list$PIP_average)) %>% select(chr)
    pos <- locus_sig_list %>% filter(PIP_average == max(locus_sig_list$PIP_average)) %>% select(pos)
    print(sentinel)
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
write.table(result,"/scratch/gen1/nnp5/Var_to_Gen_tmp/pops/results/all_result_summary_500kb_window.txt", row.names=F, quote=F, sep="\t")
write.table(result[gene_rank==1,],"/scratch/gen1/nnp5/Var_to_Gen_tmp/pops/results/prioritized_genes_500kb_window.txt", row.names=F, quote=F, sep="\t")

################# mapping features ###############
# Get gene information using biomaRt for all genes
library(biomaRt)         # Requires R/4.1.0
library(data.table)
library(tidyverse)

pops <- fread("/scratch/gen1/nnp5/Var_to_Gen_tmp/pops/results/prioritized_genes_250kb_window.txt")

ensembl <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl", version = "GRCh37")
genes <- as_tibble(getBM(attributes=c('hgnc_symbol', 'ensembl_gene_id', 'strand'),
  filters=c('ensembl_gene_id'),
  values=list(sort(unique(pops$ENSGID))),
  mart=ensembl)) %>%
  rename(gene_symbol = hgnc_symbol, gene_strand = strand)

detach("package:biomaRt", unload=TRUE)

genes <- as.data.table(genes)
setnames(genes,"ensembl_gene_id","ENSGID")
merged <- merge(pops,genes,by="ENSGID")

write.table(merged,"/scratch/gen1/nnp5/Var_to_Gen_tmp/pops/results/all_result_summary_250kb_window_gene_mapped.txt", row.names=F, quote=F, sep="\t")

##########################################################
w250 <- fread("/scratch/gen1/nnp5/Var_to_Gen_tmp/pops/results/all_result_summary_250kb_window_gene_mapped.txt") %>%
        mutate(window="250kb")
w500 <- fread("/scratch/gen1/nnp5/Var_to_Gen_tmp/pops/results/all_result_summary_500kb_window_gene_mapped.txt") %>%
        mutate(window="500kb")
w500_use <- w500[!sentinel%in%w250$sentinel,]
merged <- rbind(w250,w500_use)
merged[ENSGID=="ENSG00000182319",]$gene_symbol <- "SGK223"
merged <- merged[PoPS_score>0,]

write.table(merged,"/scratch/gen1/nnp5/Var_to_Gen_tmp/pops/results/all_results_merged_table.txt", row.names=F, quote=F, sep="\t")

merged <- fread("/scratch/gen1/nnp5/Var_to_Gen_tmp/pops/results/all_results_merged_table.txt")
signal_list <- fread("/scratch/gen1/jc824/TSH/novel_signals/TSH_signal_list.txt")
nodup <- merged[sentinel%in%signal_list$MarkerName,]
write.table(nodup,"/scratch/gen1/jc824/TSH/PoPS/all_results_merged_table.txt", row.names=F, quote=F, sep="\t")
