sink(stderr())

args <- commandArgs(T)
argc <- length(args)

if (argc < 1) {
    cat("Usage: sentinels.R HITS [ WIDTH (default", WIDTH, ")]\n")
    q()
}

tier_file = args[1]
pheno = args[2]
WIDTH <- ifelse(argc > 2, as.integer(args[3]), 500000)
pval_thr <- as.numeric(args[4])

print(paste0("total genomic window used for signal selection:",WIDTH))

suppressMessages(library(data.table))
suppressMessages(library(tidyverse))


selectSentinels <- function(data.dt) {
    regions.list <- list()
    data.dt$snpid <- paste0(data.dt$chr,"_",data.dt$bp38,"_",data.dt$a1,"_",data.dt$a2)
    while(nrow(data.dt) > 0) {
        top.dt <- data.dt[ which.min(pval) ]
        regions.list[[top.dt$snpid]] <- top.dt
        data.dt <- data.dt[ chr != top.dt$chr
            | bp38 <= top.dt$bp38 - WIDTH
            | bp38 >= top.dt$bp38 + WIDTH ]
    }
    rbindlist(regions.list)
}


tier.dt <- fread(tier_file,header=T)
tier.dt <- tier.dt %>% filter(pval <= pval_thr)
sentinels.dt <- selectSentinels(tier.dt)
sentinels.dt <- sentinels.dt %>% arrange(chr)
write.table(sentinels.dt, paste0("input/rare_variant/single_rarevar_EwWAS/",pheno,"rarevar_sentinel_variants.txt"), row.names=F, quote=F, sep="\t")