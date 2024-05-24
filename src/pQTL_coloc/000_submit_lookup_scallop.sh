#!/bin/bash

#SBATCH --job-name=scallop_pqtl_lookup
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4
#SBATCH --output=/scratch/gen1/nnp5/Var_to_Gen_tmp/logerror/%x-%j.out
#SBATCH --time=72:00:00
#SBATCH --mem=24gb
#SBATCH --account=gen1
#SBATCH --export=NONE


#Rationale: look-up in pQTL in SCALLOP pQTL - script from Chiara

PQTL_DATA="/data/gen1/reference/SCALLOP"
#SNP_LIST="/data/gen1/UKBiobank_500K/severe_asthma/Noemi_PhD/data/replsugg_valid_credset.txt"
#for 3_rs778801698:
SNP_LIST="/data/gen1/UKBiobank_500K/severe_asthma/Noemi_PhD/data/replsugg_valid_credset_chr3_noMHC.txt"
PQTL_PATH="/scratch/gen1/nnp5/Var_to_Gen_tmp/scallop_pqtl"

## List of variants with columns: ID, CHROM, POS
#cat ${SNP_LIST} | awk '{print $3"_"$4"_"$5"_"$6,$3,$4}' > ${PQTL_PATH}/snps_list.txt
#for chromsome 3_rs778801698:
cat ${SNP_LIST} | awk '$1 == "3_rs778801698_49524027_50524027" {print $3"_"$4"_"$5"_"$6,$3,$4}'  > ${PQTL_PATH}/snps_list.txt


# Set up log file
touch ${PQTL_PATH}/log_pQTL_SCALLOP_analysis
echo -ne "FILENAME\tPROTEIN\tN_LOOKUP\tN_NOMINAL\tN_SIG" > ${PQTL_PATH}/log_pQTL_SCALLOP_analysis
echo >> ${PQTL_PATH}/log_pQTL_SCALLOP_analysis


# Perform look-up using Nick's awk script
for filename in ${PQTL_DATA}/*.txt.gz
do
  ## Set file name
  name=$(basename ${filename} .txt.gz)
  ## Header for output
  zcat ${filename} | head -1 > ${PQTL_PATH}/${name}.txt
  ## Perform look-up
  /home/n/nnp5/PhD/PhD_project/Var_to_Gene/src/pQTL_coloc/000_scallop_lookup.awk \
      ${PQTL_PATH}/snps_list.txt <(zcat $filename) > ${PQTL_PATH}/${name}.txt
  ## Report look-up statistics in log file
  N_total=`cat ${PQTL_PATH}/${name}.txt | tail -n +2 -q | wc -l`
  N_nominal=`cat ${PQTL_PATH}/${name}.txt | awk '$8 < 0.05 {print $0}' | wc -l`
  N_sig=`cat ${PQTL_PATH}/${name}.txt | awk '$8 < 5E-8 {print $0}' | wc -l`
  echo -ne "${name}.txt\t${name}\t${N_total}\t${N_nominal}\t${N_sig}" >> ${PQTL_PATH}/log_pQTL_SCALLOP_analysis
  echo >> ${PQTL_PATH}/log_pQTL_SCALLOP_analysis
done