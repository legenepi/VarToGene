#!/usr/bin/env Rscript

#Rationale: data prep for rare-variant analysis in the RAP - phenotype-covariate file

library(tidyverse)
library(data.table)

bridge_file <- "/data/gen1/UKBiobank/application_88144/Bridge_eids_88144_56607.csv"
pheno_cov_file <- "/data/gen1/UKBiobank_500K/severe_asthma/Noemi_PhD/data/demo_EUR_pheno_cov_broadasthma.txt"

#match ID application 56607 to obtain phenotype-covariate file with ID application 88144:

bridge <- fread(bridge_file)
bridge$FID <- bridge$eid_88144
bridge$IID <- bridge$eid_88144
pheno_cov <- fread(pheno_cov_file) %>% select(-FID,-IID)
pheno_cov <- pheno_cov %>% rename("eid_56607" = eid)

pheno_cov_88144 <- left_join(pheno_cov,bridge,by="eid_56607") %>% relocate(IID, .before = missing) %>% relocate(FID, .before = IID)

write.table(pheno_cov_88144,"/scratch/gen1/nnp5/Var_to_Gen_tmp/rare_variant/demo_EUR_pheno_cov_broadasthma_app88144.txt",
row.names = FALSE, col.names = TRUE ,quote=FALSE, sep=" ", na = "NA")