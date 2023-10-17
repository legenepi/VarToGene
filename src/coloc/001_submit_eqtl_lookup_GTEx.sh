#!/bin/bash

#SBATCH --job-name=SA_eqtl_lookup
#SBATCH --output=/scratch/gen1/nnp5/Var_to_Gen_tmp/%x-%j.out
#SBATCH --time=0:30:0
#SBATCH --mem=16gb
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=nnp5@le.ac.uk
#SBATCH --account=gen1
#SBATCH --export=NONE

cs_all="/data/gen1/UKBiobank_500K/severe_asthma/Noemi_PhD/data/replsugg_valid_credset.txt"

#intermediate files in:
tmp_path="/scratch/gen1/nnp5/Var_to_Gen_tmp/"

module load R
Rscript ./src/coloc/preprocess_cs.R $cs_all $tmp_path
Rscript ./src/coloc/001_run_eqtl_lookup_GTEx.R $TISSUE $CREDSET


# Submit jobs
#Lets try with these, add tissues of interest later on

tissue=('Stomach' 'Small Intestine Terminal Ileum' 'Lung' 'Esophagus Muscularis' 'Esophagus Gastroesophageal Junction' 'Colon Transverse' 'Colon Sigmoid' 'Artery Tibial' 'Artery Coronary','Artery Aorta')


cs=('SA_8_81292599_C_A' 'SA_6_90963614_AT_T' 'SA_5_110401872_C_T' 'SA_2_242692858_T_C' 'SA_15_67442596_T_C' 'SA_12_56435504_C_G' 'SA_11_76296671_G_A' 'SA_9_6209697_G_A' 'SA_6_32586794_G_T' 'SA_5_131885240_C_G' 'SA_3_33042712_C_T' 'SA_2_102913642_AAAAC_A' 'SA_17_38168828_G_A' 'SA_16_27359021_C_A'  'SA_15_61068954_C_T' 'SA_12_48202941_T_C' 'SA_10_9064716_C_T')


for t in ${!tissue[*]}; do

  for c in ${!cs[@]}; do

    sbatch --export=TISSUE="${tissue[t]}",CREDSET="${cs[c]}" /scratch/gen1/atw20/pain/scripts/coloc/001_submit_eqtl_lookup_GTEx.sh

  done

done