library(tidyverse)
library(scales)

v2g_full <- read_tsv("v2g_full.txt") %>%
    mutate(gene = ifelse(gene == "GPR126", "ADGRG6", gene))

v2g_prioritised <- v2g_full %>%
    group_by(gene) %>%
    summarise(n_evidence = n_distinct(evidence)) 

v2g_evidence <- v2g_full %>%
    mutate(evidence = factor(evidence, levels = c("credible",
                                                  "pQTL",
                                                  "eQTL",
                                                  "mouse_KO",
                                                  "WES",
                                                  "rare_disease",
                                                  "PoPs",
                                                  "nearest"))) %>%
    group_by(gene) %>%
    arrange(evidence, P) %>%
    slice(1) 
   
v2g <- inner_join(v2g_prioritised, v2g_evidence) %>%
    arrange(desc(n_evidence), evidence, P)

v2g_3plus <- v2g %>%
    filter(n_evidence > 2)

cat(v2g_3plus$signal, file="v2g_3plus.snps", sep="\n")

results <- c("FEV1", "FVC", "FF", "PEF") %>%
    set_names %>%
    map_df(~paste0(., "_3plus.results") %>% read_tsv, .id="pheno")

pvals <- v2g_3plus %>%
    select(gene, signal) %>%
    left_join(results, c("signal" = "MarkerName")) %>%
    group_by(gene) %>%
    mutate(risk_allele=ifelse(sign(Zscore[pheno == "FF"]) == -1, Allele1, Allele2),
           trait=sub("FF", "FEV[1]/FVC", pheno) %>% sub("FEV1", "FEV[1]", .),
           Pval=log10(`P-value`) * sign(Zscore) * ifelse(Allele1 == risk_allele, 1, -1),
           P=cut(Pval, c(-Inf,-20,-8.3,-5,-3,3,5,8.3,20,Inf),
                 labels=c(" 10^-20", " 5%*%10^-9", " 10^-5", " 10^-3", "\"\"", "10^-3", "10^-5",
                          "5%*%10^-9", "10^-20"),
                 ordered_result = TRUE)) %>%
    ungroup %>%
    bind_rows(v2g_full %>%
                  filter(gene %in% v2g_3plus$gene) %>%
                  group_by(gene, evidence) %>%
                  slice(1) %>%
                  ungroup %>%
                  mutate(P="Yes") %>%
                  select(gene, trait=evidence, P)) %>%
    select(gene, trait, P) %>%
    pivot_wider(names_from = trait, values_from = P, values_fill = "\"\"") %>%
    pivot_longer(-gene, names_to = "trait", values_to = "P") %>%
    mutate(trait=sub("credible", "fnal_credible", trait) %>%
               factor(levels=c("fnal_credible",
                                        "eQTL",
                                        "pQTL",
                                        "WES",
                                        "mouse_KO",
                                        "rare_disease",
                                        "PoPs",
                                        "nearest",
                                        "FEV[1]/FVC", "FEV[1]", "FVC", "PEF")),
           gene=factor(gene, levels=rev(unique(gene)), ordered = T),
           P=factor(P, levels=c(" 10^-20", " 5%*%10^-9", " 10^-5", " 10^-3", "\"\"", "10^-3", "10^-5",
                                "5%*%10^-9", "10^-20", "Yes"))) %>%
    arrange(gene, trait) %>%
    group_by(trait) %>%
  #  mutate(facet=c(rep(5, ceiling(n()/5)), rep(4, ceiling(n()/5)), rep(3, ceiling(n()/5)),
   #                rep(2, ceiling(n()/5)), rep(1, floor(n()/5)))) %>%
    mutate(facet=c(rep(4, ceiling(n()/4)), rep(3, ceiling(n()/4)), 
                   rep(2, ceiling(n()/4)), rep(1, floor(n()/4)))) %>%
    ungroup %>%
    mutate(trait=factor(trait, levels=c("nearest",
                                        "PoPs",
                                        "eQTL",
                                        "rare_disease",
                                        "fnal_credible",
                                        "WES",
                                        "mouse_KO",
                                        "pQTL",
                                        "FEV[1]/FVC", "FEV[1]", "FVC", "PEF")))

prev_genes <- scan("Shrine2019_genes.txt", character())

pvals %>%
    group_by(gene) %>%
    mutate(label=paste0(gene, " (", sum(P=="Yes"), ")")) %>%
    ungroup %>%
    mutate(label=factor(label, levels=unique(label), ordered = TRUE)) %>%
    ggplot(aes(trait, label)) +
    geom_tile(aes(fill=P)) +
    scale_fill_manual(values=c("grey40", "#cc0000", "#ee5555", "#ff9999", "#ffcccc", "white",
                               "#ccccff", "#9999ff", "#5555ee", "#0000cc") %>% rev,
                      labels=label_parse()) +
    facet_wrap(~facet, scales="free", ncol=4) +
    theme(strip.background = element_blank(), strip.text = element_blank(),
        axis.text.x=element_text(angle=90, vjust=0, hjust=0, face="bold"),
        axis.title=element_blank(),
        axis.text.y=element_text(hjust=1, face="italic"),
#            face=ifelse(pvals$gene %in% prev_genes, "italic", "bold")),
        axis.ticks = element_blank(),
        legend.text.align = 0,
        legend.title = element_blank()) +
    scale_x_discrete(labels=parse(text=levels(pvals$trait)), position = "top", expand = c(0,0)) +
    scale_y_discrete(expand = c(0,0)) +
    geom_vline(xintercept = seq(0.5, 11.5, 1)) +
    geom_hline(yintercept = seq(0.5, length(unique(pvals$gene)) - 0.5, 1))
