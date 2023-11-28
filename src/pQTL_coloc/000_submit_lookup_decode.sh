#!/bin/bash

#SBATCH --job-name=decode_pqtl_lookup
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4
#SBATCH --output=/scratch/gen1/nnp5/Var_to_Gen_tmp/logerror/%x-%j.out
#SBATCH --time=72:00:00
#SBATCH --mem=70gb
#SBATCH --account=gen1
#SBATCH --export=NONE


#Rationale: look-up in pQTL in deCODE pQTL - script from Chiara
ukbb38_path="/scratch/gen1/nnp5/Var_to_Gen_tmp/ukb_pqtl"
PQTL_DATA="/data/gen1/pQTL/Ferkingstad_2021"
PQTL_PATH="/scratch/gen1/nnp5/Var_to_Gen_tmp/decode_pqtl"

#USE b38 POSITION:
##NOTE: DEPENDENCIES ON SCRIPT src/pQTL_coloc/000_preprocess_cs_b38.R TO CREATE B38 FILES:
## List of variants with columns: ID (chr_posb37_a1_a2), CHROM, POS (pos b38):
awk 'NR >1 {print $3"_"$4"_"$5"_"$6, "chr"$3, $11}' \
    ${ukbb38_path}/SA_*_b38 | grep -v "^chromosome" > ${PQTL_PATH}/snps_list.txt

# Set up log file
touch ${PQTL_PATH}/log_pQTL_decode_analysis
echo -ne "FILENAME\tPROTEIN\tN_LOOKUP\tN_NOMINAL\tN_SIG" > ${PQTL_PATH}/log_pQTL_decode_analysis
echo >> ${PQTL_PATH}/log_pQTL_decode_analysis


# Perform look-up using Nick's awk script
for filename in ${PQTL_DATA}/*.txt.gz
do
  ## Set file name
  name=$(basename ${filename} .txt.gz)
  ## Header for output
  zcat ${filename} | head -1 > ${PQTL_PATH}/${name}.txt
  ## Perform look-up
  ## Perform look-up
  /home/n/nnp5/PhD/PhD_project/Var_to_Gene/src/pQTL_coloc/000_decode_lookup.awk \
      ${PQTL_PATH}/snps_list.txt <(zcat $filename) > ${PQTL_PATH}/${name}.txt
  ## Report look-up statistics in log file
  N_total=`cat ${PQTL_PATH}/${name}.txt | tail -n +2 -q | wc -l`
  N_nominal=`cat ${PQTL_PATH}/${name}.txt | awk '$10 < 0.05 {print $0}' | wc -l`
  N_sig=`cat ${PQTL_PATH}/${name}.txt | awk '$10 < 1.8E-9 {print $0}' | wc -l`
  echo -ne "${name}.txt\t${name}\t${N_total}\t${N_nominal}\t${N_sig}" >> ${PQTL_PATH}/log_pQTL_decode_analysis
  echo >> ${PQTL_PATH}/log_pQTL_decode_analysis
done
