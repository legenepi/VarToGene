#!/usr/bin/env bash

#slurms settings

module load R
#intermediate files in:
tmp_path="/scratch/gen1/nnp5/Var_to_Gen_tmp/"




#Rationale: pipeline for Variant to Gene mapping analysis. Output from each analysis: list of genes


#1.Variant annotation
Rscript src/Variant_annotation_FAVOR.R



#2.Quantitative trait loci:
#Genomic boundaries: PIP-max causal variant +/- 1000Mb
##Divide the credible sets into separate files:
##Exclude chromosome 6 - no eQTL colocalisation for this locus.
##variables needed:
tissue=('Stomach' 'Small_Intestine_Terminal_Ileum' 'Lung' 'Esophagus_Muscularis' 'Esophagus_Gastroesophageal_Junction' 'Colon_Transverse' 'Colon_Sigmoid' 'Artery_Tibial' 'Artery_Coronary' 'Artery_Aorta')
cs=('SA_8_81292599_C_A' 'SA_6_90963614_AT_A' 'SA_5_110401872_C_T' 'SA_2_242692858_T_C' 'SA_15_67442596_T_C' 'SA_12_56435504_C_G' 'SA_11_76296671_G_A' 'SA_9_6209697_G_A' 'SA_5_131885240_C_G' 'SA_3_33042712_C_T' 'SA_2_102913642_AAAAC_A' 'SA_17_38168828_G_A' 'SA_16_27359021_C_A'  'SA_15_61068954_C_T' 'SA_12_48202941_T_C' 'SA_10_9064716_C_T')
cs_all="/data/gen1/UKBiobank_500K/severe_asthma/Noemi_PhD/data/replsugg_valid_credset.txt"

Rscript ./src/coloc/preprocess_cs.R $cs_all $tmp_path

#Obtain GWASpairs from credible set regions:
for c in ${!cs[*]}; do

  sbatch --export=CREDSET="${cs[c]}" ./src/coloc/001_submit_GWASpairs.sh

done

#Create files for GTExV8:
##'Colon Transverse' and 'Colon Sigmoid' miss from /data/gen1/ACEI/colocalisation_datasets/eQTL/GTeX/
for t in ${!tissue[*]}; do

  for c in ${!cs[@]}; do

    sbatch --export=TISSUE="${tissue[t]}",CREDSET="${cs[c]}" ./src/coloc/001_submit_eqtl_lookup_GTEx.sh

  done

done


##Get the LD matrix:
##Create the file with gtex-locus pairs:
for t in ${!tissue[*]}; do Rscript ./src/coloc/002_prepare_LDinput.R "${tissue[t]}"; done

##Get LD:
sbatch ./src/coloc/002_get_LD.sh

##to see if errors in the job: grep "Error" /scratch/gen1/nnp5/Var_to_Gen_tmp/logerror/ld-81713.out


##Run the colocalisation for GTExV8:
#run for each Tissue:
tissue='Stomach'
tissue='Small_Intestine_Terminal_Ileum'
tissue='Lung'
tissue='Esophagus_Muscularis'
tussie='Esophagus_Gastroesophageal_Junction'
tissue='Colon_Transverse'
tissue='Colon_Sigmoid'
tissue='Artery_Tibial'
tissue='Artery_Coronary'
tissue='Artery_Aorta'

for c in ${!cs[*]}; do

  N=`cat ${tmp_path}${cs[c]}_${tissue}_genes.txt | wc -l`

  sbatch --array=1-${N}%20 --export=TISSUE="${tissue}",CREDSET="${cs[c]}" ./src/coloc/003_submit_coloc_susie_GTEx.sh

 sleep 5

done