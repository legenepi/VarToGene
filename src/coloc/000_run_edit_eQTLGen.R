library(tidyverse)

freq = read_tsv("/data/gen1/LF_HRC_transethnic/eQTL/eQTLgen/2018-07-18_SNP_AF_for_AlleleB_combined_allele_counts_and_MAF_pos_added.txt.gz", 
                col_select=c("hg19_chr", "hg19_pos", "AlleleA", "AlleleB", "AlleleB_all")) %>%
  mutate(ID = paste0(hg19_chr, ":", hg19_pos, "_", pmin(AlleleA, AlleleB), "_", pmax(AlleleA, AlleleB))) %>%
  select(ID, AssessedAllele_freq = AlleleB_all)

eqtl = read_tsv("/data/gen1/reference/eqtlgen/2019-12-11-cis-eQTLsFDR-ProbeLevel-CohortInfoRemoved-BonferroniAdded.txt.gz") %>%
  mutate(ID = paste0(SNPChr, ":", SNPPos, "_", pmin(AssessedAllele, OtherAllele), "_", pmax(AssessedAllele, OtherAllele))) %>%
  relocate(ID, .before="SNP") %>%
  left_join(freq) %>%
  mutate(beta = Zscore / sqrt(2 * AssessedAllele_freq * (1 - AssessedAllele_freq) * (NrSamples + Zscore^2)), 
         se   = 1 / sqrt(2 * AssessedAllele_freq * (1 - AssessedAllele_freq) * (NrSamples + Zscore^2)))

write_tsv(x=eqtl, file="/scratch/gen1/atw20/pain/results/coloc/eqtlgen/whole_blood_cis_eqtls_withAF.txt.gz")

