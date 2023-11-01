#!/usr/bin/env bash

#slurms settings

module load R
#intermediate files in:
tmp_path="/scratch/gen1/nnp5/Var_to_Gen_tmp/"




#Rationale: pipeline for Variant to Gene mapping analysis. Output from each analysis: list of genes


#1.Variant annotation
Rscript src/Variant_annotation_FAVOR.R


#2.Quantitative trait loci:
##Exclude chromosome 6 - no eQTL colocalisation for this locus.
##Genomic boundaries: PIP-max causal variant +/- 1000Mb

##Create eQTl files in hg19 for Colon Transverse, Colon Sigmoid, Skin_Not_Sun_Exposed_Suprapubic, Skin_Sun_Exposed_Lower_leg
#From GTExV8 .parquet files and in hg38.
mkdir ${tmp_path}/liftover_gtexv8
mkdir ${tmp_path}/liftover_gtexv8/bed

dos2unix src/coloc/000A_eqtl_gtex_extraction.R src/coloc/000B_eqtl_gtex_liftover.sh \
    src/coloc/000C_eqtl_gtex_conversion.R src/coloc/000A_submit_eqtl_gtex_extraction.sh src/coloc/000C_submit_eqtl_gtex_conversion.sh
chmod o+x src/coloc/000A_eqtl_gtex_extraction.R src/coloc/000B_eqtl_gtex_liftover.sh \
    src/coloc/000C_eqtl_gtex_conversion.R src/coloc/000A_submit_eqtl_gtex_extraction.sh src/coloc/000C_submit_eqtl_gtex_conversion.sh

sbatch --array=1-22 src/coloc/000A_submit_eqtl_gtex_extraction.sh

sbatch src/coloc/000B_eqtl_gtex_liftover.sh

sbatch --array=1-22 src/coloc/000C_submit_eqtl_gtex_extraction.sh


##variables needed:
tissue=('Stomach' 'Small_Intestine_Terminal_Ileum' 'Lung' 'Esophagus_Muscularis' 'Esophagus_Gastroesophageal_Junction' 'Artery_Tibial' 'Artery_Coronary' 'Artery_Aorta' 'Colon_Transverse' 'Colon_Sigmoid' 'Skin_Sun_Exposed_Lower_leg' 'Skin_Not_Sun_Exposed_Suprapubic')
cs=('SA_8_81292599_C_A' 'SA_6_90963614_AT_A' 'SA_5_110401872_C_T' 'SA_2_242692858_T_C' 'SA_15_67442596_T_C' 'SA_12_56435504_C_G' 'SA_11_76296671_G_A' 'SA_9_6209697_G_A' 'SA_5_131885240_C_G' 'SA_3_33042712_C_T' 'SA_2_102913642_AAAAC_A' 'SA_17_38168828_G_A' 'SA_16_27359021_C_A'  'SA_15_61068954_C_T' 'SA_12_48202941_T_C' 'SA_10_9064716_C_T')
cs_all="/data/gen1/UKBiobank_500K/severe_asthma/Noemi_PhD/data/replsugg_valid_credset.txt"


##Divide the credible sets into separate files:
Rscript ./src/coloc/000_preprocess_cs.R $cs_all $tmp_path

#Obtain GWASpairs from credible set regions:
#.sh will run .R script:
for c in ${!cs[*]}; do

  sbatch --export=CREDSET="${cs[c]}" ./src/coloc/001_submit_GWASpairs.sh

done

#Create files for GTExV8:
##'Colon Transverse' and 'Colon Sigmoid' miss from /data/gen1/ACEI/colocalisation_datasets/eQTL/GTeX/
##Need to create these
#.sh will run .R script:
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
##after adding Colon and Skin tissues, to see if errors in the job: grep "Error" /scratch/gen1/nnp5/Var_to_Gen_tmp/logerror/ld-153659.out


##Run the colocalisation for GTExV8:
#.sh will run .R script:
mkdir ${tmp_path}/results
mkdir ${tmp_path}/results/gtex
#run for each Tissue:
tissue='Stomach'
tissue='Small_Intestine_Terminal_Ileum'
tissue='Lung'
tissue='Esophagus_Muscularis'
tissue='Esophagus_Gastroesophageal_Junction'
tissue='Artery_Tibial'
tissue='Artery_Coronary'
tissue='Artery_Aorta'
tissue='Colon_Transverse'
tissue='Colon_Sigmoid'
tissue='Skin_Sun_Exposed_Lower_leg'
tissue='Skin_Not_Sun_Exposed_Suprapubic'

for c in ${!cs[*]}; do

  N=`cat ${tmp_path}${cs[c]}_${tissue}_genes.txt | wc -l`

  sbatch --array=1-${N}%20 --export=TISSUE="${tissue}",CREDSET="${cs[c]}" ./src/coloc/003_submit_coloc_susie_GTEx.sh

 sleep 5

done

######QUALITY CHECKS:
#to find number of genes:
#wc -l SA_*_${tissue}.txt | sed 's/_/ /g' | sort -k 3,4 -g | awk '{print $1}'
#'Stomach' 575
#'Small_Intestine_Terminal_Ileum' 635
#'Lung' 633
#'Esophagus_Muscolaris' 564
#'Esophagus_Gastroesophageal_Junction' 576
#tissue='Artery_Tibial' 573
#tissue='Artery_Coronary' 603
#tissue='Artery_Aorta' 586
#tissue='Colon_Transverse' 610
#tissue='Colon_Sigmoid' 585
#tissue='Skin_Sun_Exposed_Lower_leg' 673
#tissue='Skin_Not_Sun_Exposed_Suprapubic' 670

##Check that all genes for each tissue have been analysed:
grep ${tissue} ${tmp_path}/logerror/coloc_susie_gtex*.out | awk -F ":" '{print $1}' | sort -u | wc -l

##Check how many genes per tissue have been analysed for colocalisation:
ls -lthr  ${tmp_path}/results/gtex/*all_coloc.rds | grep ${tissue} | wc -l
ls -lthr ${tmp_path}/results/gtex/*all_susie*.rds | grep ${tissue} | wc -l

#Find statistically significant colocalisation results and add results into var2gene_raw.xlsx:
Rscript ./src/coloc/004_concat_coloc_results.R