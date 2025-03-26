#!/bin/bash

#SBATCH --job-name=decode_pqtl_lookup_additionalcredsetMarch2025
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4
#SBATCH --output=/scratch/gen1/nnp5/Var_to_Gen_tmp/logerror/%x-%j.out
#SBATCH --time=72:00:00
#SBATCH --mem=70gb
#SBATCH --account=gen1
#SBATCH --export=NONE

#call this script as:
#From Chiara/Kayesha script and Nick:
sbatch src/pQTL_coloc/000_submit_lookup_decode.sh

#Rationale: look-up in pQTL in deCODE pQTL - script from Chiara
PQTL_DATA="/data/gen1/pQTL/Ferkingstad_2021"
#mkdir /scratch/gen1/nnp5/Var_to_Gen_tmp/decode_pqtl
PQTL_PATH="/scratch/gen1/nnp5/Var_to_Gen_tmp/decode_pqtl"

#USE b38 POSITION:
#from the Credset.xlsx file create a file with:
#ID (chr_posb37_a1_a2)	Chr_N	Pos_b38
#2_102926362_A_G	chr2	102309902
#cp /home/n/nnp5/PhD/PhD_project/Var_to_Gene/input/Additional_credset_snps_March2025/Credset_snps_decode_input.txt ${PQTL_PATH}/snps_list.txt

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

##filter out gene name with significant pQTL:
#awk 'NR > 1 && $5 > 0 {print $1}' ${PQTL_PATH}/log_pQTL_decode_analysis | sed 's/.txt//g' \
#    > /home/n/nnp5/PhD/PhD_project/Var_to_Gene/input/Additional_credset_snps_March2025/decode_pqtl_var2genes_raw_additional_credsetMarch2025
##no results for decode