#!/usr/bin/env Rscript

#Rationale: Find recombination boundaries in which replicated sentinel variants are mapped.
#Recombination boundaries will be used as boundaries for the colocalisation analysis.

rm(list= ls())

library(data.table)
library(tidyverse)

args = commandArgs(trailingOnly = TRUE)
recom_bnd_file = args[1] #/home/n/nnp5/PhD/PhD_project/Var_to_Gene/input/boundaries_recomb_hotspots_ldetect.txt
sugg_sent_file = args[2] #/data/gen1/UKBiobank_500K/severe_asthma/Noemi_PhD/data/fine_mapping_regions_replicated_suggestive_input

recom_bnd <- fread(args[1])
sugg_bnd <- fread(args[2])

recom_bnd <- fread("/home/n/nnp5/PhD/PhD_project/Var_to_Gene/input/boundaries_recomb_hotspots_ldetect.txt")
colnames(recom_bnd) <- c("chr","start_bnd","stop_bnd")
sugg_sent <- fread("/data/gen1/UKBiobank_500K/severe_asthma/Noemi_PhD/data/fine_mapping_regions_replicated_suggestive_input")
colnames(sugg_sent) <- c('snpid','chr','bp','start','end')
#put together sentinel on chromosome 17, it went splitted because chr:pos_A1_A2:
sugg_sent <-  sugg_sent %>% separate(snpid, sep = "_", into = c("snpid_1", "snpid_2", "snpid_3"), extra = "drop")
sugg_sent <-  sugg_sent %>% separate(bp, sep = "_", into = c("bp_1", "bp_2", "bp_3"), extra = "drop")
sugg_sent[18,1] <- paste0(sugg_sent$snpid_1[18],"_",sugg_sent$snpid_2[18],"_",sugg_sent$snpid_3[18])
sugg_sent$chr <- paste0("chr",sugg_sent$chr)

#divide loci with single or multiple sentinel vars:
sugg_sent_single <- sugg_sent %>% filter(is.na(snpid_2) | snpid_2 == "CCG")
sugg_sent_mult <- sugg_sent %>% filter(!is.na(snpid_2) & snpid_2 != "CCG")

#initialise the dataframe with: snpid_1  chr      bp_1     start       end start_bnd  stop_bnd
columns <- c("snpid_1","chr","bp_1","start","end","start_bnd","stop_bnd")
snp_recom_bnd <- data.frame(matrix(nrow = 0, ncol = length(columns)))
colnames(snp_recom_bnd) <- columns

for (i in 1:dim(sugg_sent_single)[1]){
    snpid_tmp <- sugg_sent_single$snpid_1[i]
    chr_tmp <- sugg_sent_single$chr[i]
    bp_tmp <- as.numeric(sugg_sent_single$bp_1[i])
    recom_bnd_tmp <- recom_bnd %>% filter(chr == chr_tmp & start_bnd <= bp_tmp & bp_tmp <= stop_bnd)
    sentinel <- sugg_sent_single %>% filter(snpid_1 == snpid_tmp)
    snp_recom_bnd_tmp <- left_join(sentinel, recom_bnd_tmp, by="chr") %>% select(!c(snpid_2,snpid_3,bp_2,bp_3))
    snp_recom_bnd <- rbind(snp_recom_bnd, snp_recom_bnd_tmp)
    }

#TO BE FINISHED:
for (i in 1:dim(sugg_sent_mult)[1]){
    snpid_tmp <- sugg_sent_single$snpid_1[i]
    chr_tmp <- sugg_sent_single$chr[i]
    bp_tmp <- as.numeric(sugg_sent_single$bp_1[i])
    recom_bnd_tmp <- recom_bnd %>% filter(chr == chr_tmp & start_bnd <= bp_tmp & bp_tmp <= stop_bnd)
    sentinel <- sugg_sent_single %>% filter(snpid_1 == snpid_tmp)
    snp_recom_bnd_tmp <- left_join(sentinel, recom_bnd_tmp, by="chr")
    ##recreate snpid with both snpid as well as bp:

    snp_recom_bnd <- rbind(snp_recom_bnd, snp_recom_bnd_tmp)
    }





