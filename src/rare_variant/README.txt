#Rationale: Micellaneous knowldege about dx toolkit, RAP and CLI to run the rare-variant single and gene-based collapsing analysis analysis step by step

#In REGENIE, phenotype and covariates column names:
file: /rfs/TobinGroup/data/UKBiobank/application_88144/demo_EUR_pheno_cov_broadasthma_app88144.txt
pheno="broad_pheno_1_5_ratio"
--covarColList age_at_recruitment,age2,PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10,genetic_sex

#Download the CLI:
https://documentation.dnanexus.com/downloads

#Tutorial videos at:
https://dnanexus.gitbook.io/uk-biobank-rap/working-on-the-research-analysis-platform/accessing-and-using-the-command-line-interface

#When indexing files on the RAP in the CLI: always put the project ID, it makes it run faster:
project-GGzFY70JBJzVx22v4Yj980J1:/Bulk/Genotype Results/CEL files/11/

#Upload the phenotype file on the platform through dx toolkit:
dx login
dx logout
dx upload /rfs/TobinGroup/data/UKBiobank/application_88144/demo_EUR_pheno_cov_broadasthma_app88144.txt
dx describe #description of the tool