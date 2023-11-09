sink(stderr())
    
args <- commandArgs(T)
argc <- length(args)

batch <- as.integer(args[1])

suppressMessages(library(data.table))
suppressMessages(library(tidyverse))
suppressMessages(library(dplyr))
suppressMessages(library(reticulate))

list <- fread("/scratch/gen1/jc824/transethnic/eQTL/lung_eQTL/lung_eQTL_lookup_result.txt",header=T)
n <- nrow(list)
i1=batch*10000+1
i2=batch*10000+10000
if(i2>n){
    i2=n
}
#gene_probe <- fread("/data/gen1/reference/lung_eQTL/tabMerged_anno.txt",header=T)
#gene_probe$ProbeSet <- gsub("_at","",gene_probe$ProbeSet)
#setnames(gene_probe,"gene","gene_id")
eQTL_lookup <- data.frame(matrix(ncol = 5,nrow = 0))

for(i in i1:i2){
   print(i)
   sentinel <- list$sentinel[i]
   chr <- list$chr[i]
   pos <- list$pos[i]
   trait <- list$trait[i]
   probe <- list$ProbeSet[i] 
   start <- pos-1000000
   if(start<0){
       start=0
   }
   end <- pos+1000000
   gwas <- fread(paste0("/scratch/gen1/jc824/transethnic/conditional_LF_signals_EUR/",trait,"_",sentinel,"_conditional.txt",sep=""),header=T)
   gwas <- gwas[gwas$BP>=start,]
   gwas <- gwas[gwas$BP<=end,]
   eQTL <- fread(paste0("/scratch/gen1/jc824/transethnic/eQTL/lung_eQTL/eQTL_region_stat/",sentinel,"_eQTL_region.txt",sep=""),header=T)
   setnames(eQTL,"#Probe","ProbeSet")
   eQTL <- eQTL[eQTL$ProbeSet==probe,]
   eQTL <- mutate(eQTL,SNP=paste(CHR,BP,pmin(toupper(Allele1),toupper(Allele2)),pmax(toupper(Allele1),toupper(Allele2)),sep="_"))
   overlapped <- merge(gwas,eQTL,by="SNP",all.x=F,all.y=F)
   if(pos>=min(overlapped$BP.x)&(pos<=max(overlapped$BP.x))){
        bim <- fread(paste0("/scratch/gen1/jc824/transethnic/eQTL/LD_matrix/",sentinel,"_1mb.bim"),header=F)
        setnames(bim,"V2","SNP")
        merged <- merge(overlapped,bim,by="SNP",all.x=F,all.y=F)
        if(pos>=min(merged$BP.x)&(pos<=max(merged$BP.x))){
            eQTL_lookup <<- rbind(eQTL_lookup,list[i,])
            merged <- merged[order(BP.x)]
            snp.use <- select(merged,"SNP")
            gwas.stat <- select(merged,"SNP","freq","b","se","p","N")
            gwas.stat$freq <- as.numeric(gwas.stat$freq)
            gwas.stat$b <- as.numeric(gwas.stat$b)
            gwas.stat$se <- as.numeric(gwas.stat$se)
            gwas.stat$p <- as.numeric(gwas.stat$p)
            gwas.stat$N <- as.numeric(gwas.stat$N)
            gwas.stat <- mutate(gwas.stat,maf=ifelse(freq>0.5,1-freq,freq))
            gwas.stat <- mutate(gwas.stat,beta=ifelse(freq>0.5,-1*b,b))
            gwas.stat <- select(gwas.stat,"SNP","maf","beta","se","p","N")
            setnames(gwas.stat,"beta","b")
            eQTL.stat <- select(merged,"SNP","Freq1","Effect","StdErr","P")
            setnames(eQTL.stat,"Freq1","freq")
            setnames(eQTL.stat,"Effect","b")
            setnames(eQTL.stat,"StdErr","se")
            setnames(eQTL.stat,"P","p")
            eQTL.stat$freq <- as.numeric(eQTL.stat$freq)
            eQTL.stat$b <- as.numeric(eQTL.stat$b)
            eQTL.stat$se <- as.numeric(eQTL.stat$se)
            eQTL.stat$z <- eQTL.stat$b/eQTL.stat$se
            eQTL.stat$N <- 1/(2*eQTL.stat$freq*(1-eQTL.stat$freq)*eQTL.stat$se*eQTL.stat$se)-eQTL.stat$z*eQTL.stat$z
            eQTL.stat <- mutate(eQTL.stat,maf=ifelse(freq>0.5,1-freq,freq))
            eQTL.stat <- mutate(eQTL.stat,beta=ifelse(freq>0.5,-1*b,b))         
            eQTL.stat <- select(eQTL.stat,"SNP","maf","beta","se","p","N")
            setnames(eQTL.stat,"beta","b")
            write.table(snp.use,file=paste0("/scratch/gen1/jc824/transethnic/eQTL/lung_eQTL/coloc_input/",sentinel,"_",trait,"_",probe,"_snp.use.txt",sep=""),col.names=F, row.names = F, sep = "\t", quote=F)
            write.table(gwas.stat,file=paste0("/scratch/gen1/jc824/transethnic/eQTL/lung_eQTL/coloc_input/",sentinel,"_",trait,"_",probe,"_gwas.stat.txt",sep=""),col.names=T, row.names = F, sep = "\t", quote=F)
            write.table(eQTL.stat,file=paste0("/scratch/gen1/jc824/transethnic/eQTL/lung_eQTL/coloc_input/",sentinel,"_",trait,"_",probe,"_eQTL.stat.txt",sep=""),col.names=T, row.names = F, sep = "\t", quote=F) 
        }
    }
}

sink()
write.table(eQTL_lookup, "", col.names=T, row.names = F, sep = "\t", quote=F)
