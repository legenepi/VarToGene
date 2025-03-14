#!/usr/bin/env bash

#Rationale: run PoPS with severe asthma summary statistics
#Based on scripts and log files found in: /data/gen1/TSH/PoPS
#PoPS v0.1: https://github.com/FinucaneLab/pops/tree/add-license-1
#software stored in /data/gen1/LF_HRC_transethnic/PoPS/pops

module load R
#intermediate files in:
tmp_path="/scratch/gen1/nnp5/Var_to_Gen_tmp/"

#Download MAGMA and PoPS:
cd /home/n/nnp5/software
wget https://ctg.cncr.nl/software/MAGMA/prog/magma_v1.10.zip
unzip https://ctg.cncr.nl/software/MAGMA/prog/magma_v1.10.zip
cd /home/n/nnp5/PhD/PhD_project/Var_to_Gene

#Step 0: Generate gene-level association statistics as by MAGMA to have genes.out and genes.raw files:
##bfile and gene-annot: already downloaded with respect to TSH/lung function analysis - use the same.
##pval: sumstats file from severe asthma GWAS as this:
#SNP	p	N
#rs147324274	0.9139	124358
#rs369318156	0.4323	124358

#Create the .sumstats file:
awk '{print $1, $9, $3=46086}' \
    /data/gen1/UKBiobank_500K/severe_asthma/Noemi_PhD/data/maf001_broad_pheno_1_5_ratio_betase_input_mungestat | \
    sed 's/snpid pval 46086/SNP p N/g' > ${tmp_path}/pops/SA.sumstat

/home/n/nnp5/software/magma \
    --bfile /data/gen1/LF_HRC_transethnic/PoPS/data/1000G.EUR \
	  --gene-annot /data/gen1/LF_HRC_transethnic/PoPS/data/magma_0kb.genes.annot \
	  --pval ${tmp_path}/pops/SA.sumstat ncol=N \
	  --gene-model snp-wise=mean \
	  --out ${tmp_path}/pops/SA


#Step 1: Select features
#the features are in: /data/gen1/LF_HRC_transethnic/PoPS/data/PoPS.features.txt.gz
sbatch src/PoPS/submit_pops.sh


