#!/bin/bash

#Rationale: run Single-variant rare-variant ExWAS in the RAP from command-line interface CLI, using dx toolkit syntax
#Command-line interface (script)
#Script form Kath:
#What follows is my bash script for running a GWAS (using regenie) on the command-line interface.  I used the genotyping data for the first step and the exome data for the second.  For more information on the command-line interface and how to use it, I'd recommend the videos on this page: https://dnanexus.gitbook.io/uk-biobank-rap/working-on-the-research-analysis-platform/accessing-and-using-the-command-line-interface
#https://dnanexus.gitbook.io/uk-biobank-rap/getting-started/research-analysis-platform-training-webinars

cd /home/n/nnp5/PhD/PhD_project/Var_to_Gene/

#Start of script for regenie GWAS:
dx login
# You will need to enter your username and password at this stage, and select the number corresponding to your project from the list
dx select Severe_asthma
dx mkdir analysis
data_file_dir="/analysis/"

# Merging chromosomes together (genotyping data):
#The RAP folders names have space on that ! so need to use '\' to specify that:
# demo_EUR_pheno_cov_broadasthma_app88144.txt is the phenotype file I generated on ALICE3.
# running locally in the CLI, and then it gives a job that I am able to monitor on the RAP.
run_merge="cp /mnt/project/Bulk/Genotype\ Results/Genotype\ calls/ukb22418_c[1-9]* . ; ls *.bed | sed -e 's/.bed//g'> files_to_merge.txt; plink --merge-list files_to_merge.txt --make-bed --autosome-xy --out ukb22418_c1_22_v2_merged; rm files_to_merge.txt;"
dx run swiss-army-knife -iin="/demo_EUR_pheno_cov_broadasthma_app88144.txt" -icmd="${run_merge}" --tag="Merge" --instance-type "mem1_ssd1_v2_x16" --destination="${data_file_dir}" --brief --yes

# QC-ing the merged genotyping data:
awk {'print $1, $2'} /rfs/TobinGroup/data/UKBiobank/application_88144/demo_EUR_pheno_cov_broadasthma_app88144.txt | tail -n +2 > /scratch/gen1/nnp5/Var_to_Gen_tmp/rare_variant/ukb_eur_ids
dx cd analysis
dx upload ukb_eur_ids
run_plink_qc="plink2 --bfile ukb22418_c1_22_v2_merged --keep ukb_eur_ids --autosome --maf 0.01 --mac 20 --geno 0.1 --hwe 1e-15 --mind 0.1 --write-snplist --write-samples --no-id-header --out snp_qc_pass"
dx run swiss-army-knife -iin="${data_file_dir}/ukb22418_c1_22_v2_merged.bed" -iin="${data_file_dir}/ukb22418_c1_22_v2_merged.bim" -iin="${data_file_dir}/ukb22418_c1_22_v2_merged.fam" -iin="${data_file_dir}/ukb_eur_ids" -icmd="${run_plink_qc}" --tag="plink_QC_eur_ukb" --instance-type "mem1_ssd1_v2_x16" --destination="${data_file_dir}" --brief --yes

#  Running regenie step 1:
##Using data that represents the whole genome, either imputed or whole genome seq:
run_regenie_step1="regenie --step 1 \
--lowmem \
--out sa_results \
--bed ukb22418_c1_22_v2_merged \
--phenoFile demo_EUR_pheno_cov_broadasthma_app88144.txt \
--covarFile demo_EUR_pheno_cov_broadasthma_app88144.txt \
--extract snp_qc_pass.snplist \
--phenoCol broad_pheno_1_5_ratio \
--covarCol age_at_recruitment,age2,PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10,genetic_sex \
--bsize 1000 \
--bt \
--loocv \
--gz \
--threads 16"
dx run swiss-army-knife -iin="${data_file_dir}/ukb22418_c1_22_v2_merged.bed" -iin="${data_file_dir}/ukb22418_c1_22_v2_merged.bim" -iin="${data_file_dir}/ukb22418_c1_22_v2_merged.fam" -iin="/demo_EUR_pheno_cov_broadasthma_app88144.txt" -iin="${data_file_dir}/snp_qc_pass.snplist" -icmd="${run_regenie_step1}" --tag="Step1" --instance-type "mem1_ssd1_v2_x16" --destination="${data_file_dir}" --brief --yes

# QC-ing the exome data:
exome_file_dir="/Bulk/Exome sequences/Population level exome OQFE variants, PLINK format - final release"
field_name="ukb23158"
for chr in {1..22}
do
  run_plink_wes="plink2 \
  --bfile ukb23158_c${chr}_b0_v1 \
  --no-psam-pheno \
  --keep ukb_eur_ids \
  --autosome \
  --mac 3 \
  --geno 0.1 \
  --hwe 1e-15 \
  --mind 0.1 \
  --write-snplist \
  --write-samples \
  --no-id-header \
  --out WES_c${chr}_snps_qc_pass"
  dx run swiss-army-knife \
  -iin="${exome_file_dir}/${field_name}_c${chr}_b0_v1.bed" \
  -iin="${exome_file_dir}/${field_name}_c${chr}_b0_v1.bim" \
  -iin="${exome_file_dir}/${field_name}_c${chr}_b0_v1.fam" \
  -iin="${data_file_dir}/ukb_eur_ids" \
  -icmd="${run_plink_wes}" \
  --tag="QC_ExWAS" \
  --instance-type "mem1_ssd1_v2_x16" \
  --name "QC_chr${chr}" \
  --destination="${data_file_dir}" \
  --brief \
  --yes
done

# Running regenie step 2:
exome_file_dir="/Bulk/Exome sequences/Population level exome OQFE variants, PLINK format - final release"
EXCLUDE="ukb23158_500k_OQFE.90pct10dp_qc_variants.txt"
field_name="ukb23158"
for chr in {1..22}
do
  run_regenie_cmd="regenie \
  --step 2 \
  --bed ukb23158_c${chr}_b0_v1 \
  --out ExWAS_SA_assoc.c${chr} \
  --phenoFile demo_EUR_pheno_cov_broadasthma_app88144.txt \
  --covarFile demo_EUR_pheno_cov_broadasthma_app88144.txt \
  --bt \
  --approx \
  --firth-se \
  --firth \
  --extract WES_c${chr}_snps_qc_pass.snplist \
  --exclude $EXCLUDE \
  --phenoCol broad_pheno_1_5_ratio \
  --covarCol age_at_recruitment,age2,PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10,genetic_sex \
  --bsize 400 \
  --pred sa_results_pred.list \
  --pThresh 0.05 \
  --minMAC 3 \
  --threads 16 \
  --gz"
  dx run swiss-army-knife \
  -iin="${exome_file_dir}/${field_name}_c${chr}_b0_v1.bed" \
  -iin="${exome_file_dir}/${field_name}_c${chr}_b0_v1.bim" \
  -iin="${exome_file_dir}/${field_name}_c${chr}_b0_v1.fam" \
  -iin="${data_file_dir}/WES_c${chr}_snps_qc_pass.snplist" \
  -iin="${exome_file_dir}/helper_files/${EXCLUDE}" \
  -iin="/demo_EUR_pheno_cov_broadasthma_app88144.txt" \
  -iin="${data_file_dir}/sa_results_pred.list" \
  -iin="${data_file_dir}/sa_results_1.loco.gz" \
  -icmd="${run_regenie_cmd}" \
  --tag="Step2" \
  --instance-type "mem1_ssd1_v2_x16" \
  --name="regenie_step2_chr${chr}" \
  --destination="${data_file_dir}" \
  --brief \
  --yes
done

#Downlaod a log file to have it locally on ALICE3:
#Var_to_Gene/input/regenie_step2_chr5:
#REGENIE v3.1.1.gz
#* case-control counts for each trait:
#'broad_pheno_1_5_ratio': 7413 cases and 36955 controls

# Merging output of regenie and formatting for visualising using LocusZoom:
merge_cmd='out_file="assoc.regenie.merged.txt"; cp /mnt/project/analysis/*.regenie.gz . ; gunzip *.regenie.gz ; echo -e "#CHROM\tGENPOS\tID\tALLELE0\tALLELE1\tA1FREQ\tN\tTEST\tBETA\tSE\tCHISQ\tLOG10P\tEXTRA" > $out_file ; files="./*.regenie"; for f in $files; do tail -n+2 $f | tr " " "\t" >> $out_file; done ; rm *.regenie'
dx run swiss-army-knife \
    -iin="${data_file_dir}/ExWAS_SA_assoc.c1_broad_pheno_1_5_ratio.regenie.gz" \
    -icmd="${merge_cmd}" \
    --tag="Merge_regenie_results" \
    --instance-type "mem1_ssd1_v2_x16" \
    --destination="${data_file_dir}" \
    --brief \
    --yes

#Run LocusZoom:
dx run locuszoom -igwas_files="${data_file_dir}assoc.regenie.merged.txt"

# End of script.

