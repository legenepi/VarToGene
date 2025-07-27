#!/usr/bin/env bash

#Rationale: do eQTL colocalisation using the final list of signal variants and credible set - July 2025.

module load R
#intermediate files in:
tmp_path="/scratch/gen1/nnp5/Var_to_Gen_tmp/"

################
#2.1 COLOCALISATION - eQTL
###############
#GTExV8
################

##copy the liftOver files for from Kath's scratch for the chromosomes, for tissues:
#"Colon_Transverse", "Colon_Sigmoid", "Skin_Sun_Exposed_Lower_leg", "Skin_Not_Sun_Exposed_Suprapubic":
for i in {2,5,6,8,9,10,11,12,15}; do
cp /scratch/ukb/kaf19/Noemi_V2G/eqtl_gtex/processing/liftover_gtexv8/*v8.EUR.allpairs.chr$i.hg19.txt.gz \
    $tmp_path/liftover_gtexv8/
done

##variables needed:
cs=('SA_2_102926362_G_A' 'SA_2_242692858_C_T' 'SA_5_110401872_T_C' 'SA_5_110404999_A_C' 'SA_5_131887986_C_A'
'SA_5_131819921_C_A' 'SA_6_90963614_A_AT' 'SA_8_81292599_A_C' 'SA_9_6209697_A_G' 'SA_10_9049253_C_T' 'SA_10_9064361_T_C'
'SA_11_76296671_A_G' 'SA_12_56435504_G_C' 'SA_12_57493727_T_G' 'SA_15_67442596_C_T')

#Manually create the below file, with PIPs:
#nano input/Additional_credset_snps_March2025/credset_July2025.txt
#the colnames to match in the ./src/coloc/000_preprocess_cs_chr12_chr5_March2025.R
cs_all="input/Additional_credset_snps_March2025/credset_July2025.txt"

##Divide the credible sets into separate files:
Rscript ./src/coloc/000_preprocess_cs_chr12_chr5_March2025.R $cs_all $tmp_path

#Obtain GWASpairs from credible set regions:
#.sh will run .R script:
for c in ${!cs[*]}; do

  sbatch --export=CREDSET="${cs[c]}" ./src/coloc/001_submit_GWASpairs.sh

done

#Create files for GTExV8:
##'Colon Transverse' and 'Colon Sigmoid' miss from /data/gen1/ACEI/colocalisation_datasets/eQTL/GTeX/
##Need to create these
#.sh will run .R script:
tissue=('Stomach' 'Small_Intestine_Terminal_Ileum' 'Lung' 'Esophagus_Muscularis' 'Esophagus_Gastroesophageal_Junction'
'Artery_Tibial' 'Artery_Coronary' 'Artery_Aorta' 'Colon_Transverse' 'Colon_Sigmoid' 'Skin_Sun_Exposed_Lower_leg'
'Skin_Not_Sun_Exposed_Suprapubic')
#need to rename the cs as they are saved with a different allele order:
cs=('SA_2_102926362_G_A' 'SA_2_242692858_C_T' 'SA_5_110401872_T_C' 'SA_5_110404999_A_C' 'SA_5_131887986_C_A'
'SA_5_131819921_C_A' 'SA_6_90963614_A_AT' 'SA_8_81292599_A_C' 'SA_9_6209697_A_G' 'SA_10_9049253_C_T' 'SA_10_9064361_T_C'
'SA_11_76296671_A_G' 'SA_12_56435504_G_C' 'SA_12_57493727_T_G' 'SA_15_67442596_C_T')

for t in ${!tissue[*]}; do

  for c in ${!cs[@]}; do

    sbatch --export=TISSUE="${tissue[t]}",CREDSET="${cs[c]}" ./src/coloc/001_submit_eqtl_lookup_GTEx.sh

  done

done

#To see if error: grep "Error" logerror/SA_eqtl_lookup-*
#they should be 180 each (15 credset X 12 tissues):
#ls -lthr | grep "genes.txt" | wc -l
#ls -lthr | grep "pairs.txt.gz" | wc -
#ls -lthr | grep "lookup.txt.gz" | wc -l

#Save the look-up with both GWAS and eQTL p-value <= 5 x 10-6:
zcat $tmp_path/SA_*_lookup.txt.gz | head -n1 > $tmp_path/header_gtex_lookup
zcat $tmp_path/SA_*_lookup.txt.gz | awk '$13 <= 0.000005 && $24 <= 0.000005' \
    > $tmp_path/gtex_lookup_suggestive_july2025.txt

cat $tmp_path/header_gtex_lookup $tmp_path/gtex_lookup_suggestive_july2025.txt > input/gtex_lookup_suggestive_july2025.txt
cp input/gtex_lookup_suggestive_july2025.txt /rfs/TobinGroup/nnp5/data/

##Get the LD matrix:
##Create the file with gtex-locus pairs:
for t in ${!tissue[*]}; do Rscript ./src/coloc/002_prepare_LDinput.R "${tissue[t]}"; done

##Get LD:
#with parameters for GTExV8
sbatch ./src/coloc/002_get_LD.sh

##to see if errors in the job: grep "Error" /scratch/gen1/nnp5/Var_to_Gen_tmp/logerror/ld-4744247.out

##Run the colocalisation for GTExV8:
#.sh will run .R script:
mkdir ${tmp_path}/results
mkdir ${tmp_path}/results/gtex

 TO BE FINISHED
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

  N=`cat ${tmp_path}/${cs[c]}_${tissue}_genes.txt | wc -l`

  sbatch --array=1-${N}%20 --export=TISSUE="${tissue}",CREDSET="${cs[c]}" ./src/coloc/003_submit_coloc_susie_GTEx.sh

 sleep 5

done

######QUALITY CHECKS:
##Check that all genes for each tissue have been analysed:
for c in ${!tissue[*]}; do
echo ${tissue[c]}; echo "Total:"
wc -l ${tmp_path}/SA_*_${tissue[c]}_genes.txt | sed 's/_/ /g' | sort -k 3,4 -g | awk '{print $1}'
echo "Analysed:"; grep ${tissue[c]} ${tmp_path}/logerror/coloc_susie_gtex*.out | awk -F ":" '{print $1}' | sort -u | wc -l
##Check how many genes per tissue have been analysed for colocalisation:
echo "Coloc analysed genes:"
ls -lthr  ${tmp_path}/results/gtex/*all_coloc.rds | grep ${tissue} | wc -l
echo "Coloc susie analysed genes:"
ls -lthr ${tmp_path}/results/gtex/*all_susie*.rds | grep ${tissue} | wc -l
done


######################################################
###eqtlGen eQTL###
######################################################
mkdir ${tmp_path}/results/eqtlgen
mkdir ${tmp_path}/eqtlgen

sbatch src/coloc/000_submit_edit_eQTLGen.sh

#Create files for eQTLGen colocalisation:
for c in ${!cs[@]}; do

    sbatch --export=CREDSET="${cs[c]}" ./src/coloc/001_submit_eqtl_lookup_eQTLGen.sh

done

#save lookup gwas and eqtl p-value <= 0.000005:
zcat $tmp_path/eqtlgen/SA_*_lookup.txt.gz | head -n1 > $tmp_path/header_eqtlgen_lookup
zcat $tmp_path/eqtlgen/SA_*_eqtlGenWB_lookup.txt.gz  | awk '$13 <= 0.000005 && $15 <= 0.000005' | \
    grep -v "Replicated" \
    > $tmp_path/eqtlgen_lookup_suggestive_july2025.txt

cat $tmp_path/header_eqtlgen_lookup $tmp_path/eqtlgen_lookup_suggestive_july2025.txt > input/eqtlgen_lookup_suggestive_july2025.txt
cp input/eqtlgen_lookup_suggestive_july2025.txt /rfs/TobinGroup/nnp5/data/

##Get the LD matrix:
##Create the file with gtex-locus pairs:
Rscript src/coloc/002_prepare_LDinput_eqtlgen.R "eqtlGenWB"

##Get LD:
#with parameters for eqtlGen
sbatch src/coloc/002_get_LD.sh
##to see if errors in the job: grep "Error" /scratch/gen1/nnp5/Var_to_Gen_tmp/logerror/ld-4745617.out

##Run the colocalisation for eQTlGen:
#.sh will run .R script:
#need to rename the cs as they are saved with a different allele order:
cs=('SA_2_102926362_G_A' 'SA_2_242692858_C_T' 'SA_5_110401872_T_C' 'SA_5_110404999_A_C' 'SA_5_131887986_C_A'
'SA_5_131819921_C_A' 'SA_6_90963614_A_AT' 'SA_8_81292599_A_C' 'SA_9_6209697_A_G' 'SA_10_9049253_C_T' 'SA_10_9064361_T_C'
'SA_11_76296671_A_G' 'SA_12_56435504_G_C' 'SA_12_57493727_T_G' 'SA_15_67442596_C_T')

for c in ${!cs[*]}; do

  N=`cat ${tmp_path}eqtlgen/${cs[c]}_eqtlGenWB_genes.txt | wc -l`

  sbatch --array=1-${N}%20 --export=CREDSET="${cs[c]}" ./src/coloc/003_submit_coloc_susie_eQTLGen.sh

 sleep 5

done

######QUALITY CHECKS:
#to find number of genes: 468
wc -l ${tmp_path}/eqtlgen/SA_*_eqtlGenWB_genes.txt | sed 's/_/ /g' | sort -k 3,4 -g | awk '{print $1}'
##Check that all genes for each tissue have been analysed:
grep "eqtlGenWB" ${tmp_path}/logerror/coloc_susie_eqtlgen*.out | awk -F ":" '{print $1}' | sort -u | wc -l

##Check how many genes per tissue have been analysed for colocalisation:
ls -lthr  ${tmp_path}/results/eqtlgen/*all_coloc.rds | grep "eqtlGenWB" | wc -l
ls -lthr ${tmp_path}/results/eqtlgen/*all_susie*.rds | grep "eqtlGenWB" | wc -l

#######

#gtex converted in gene symbol:
#https://www.biotools.fr/human/ensembl_symbol_converter
#added in GTExV8_eQTL_genes_symbol table in the input/var2gene.xlsx file.


#Find statistically significant colocalisation results for GTExV8 and eqtlGen eQTL, and add results into var2gene_raw.xlsx:
Rscript ./src/coloc/004_concat_coloc_results_additionalcredset_March2025.R


#################################################
###UBCLung eQTL###
#################################################
mkdir ${tmp_path}/results/ubclung
mkdir ${tmp_path}/ubclung

for c in ${!cs[*]}; do

  sbatch --export=CREDSET="${cs[c]}" ./src/coloc_UBClung/000_submit_lookup_lung_eQTL.sh

done

##Get the LD matrix:
##Create the file with gtex-locus pairs:
Rscript src/coloc_UBClung/002_prepare_LDinput_ubclung.R "UBCLung"

##Get LD:
#with parameters for UBCLung
sbatch src/coloc/002_get_LD.sh
##to see if errors in the job: grep "Error" /scratch/gen1/nnp5/Var_to_Gen_tmp/logerror/ld-4758076.out

##Create UBCLung eQTL regional data from eQTL summary stats:
mkdir ${tmp_path}/ubclung/eQTL_region_stat/
module load tabix
lung_eQTL="/data/gen1/reference/lung_eQTL"
touch ${tmp_path}/ubclung/ubclung_lookup_suggestive_july2025.txt

for c in ${!cs[@]}
    do
      CREDSET="${cs[c]}"
      read chr < <(echo $CREDSET | awk -F "_" '{print $2}')
      read pos < <(echo $CREDSET | awk -F "_" '{print $3}')
      pos1=`expr $pos - 500000`
      pos2=`expr $pos + 500000`
      if (( pos1 < 0 )); then
          pos1=0
      fi
      tabix -h ${lung_eQTL}/METAANALYSIS_Laval_UBC_Groningen_chr${chr}_formatted.txt.gz ${chr}:${pos1}-${pos2} > ${tmp_path}/ubclung/eQTL_region_stat/${CREDSET}_eQTL_region.txt
      echo $CREDSET >> ${tmp_path}/ubclung/ubclung_lookup_suggestive_july2025.txt
      awk '$10 < 0.000005 {print}' ${tmp_path}/ubclung/eQTL_region_stat/${CREDSET}_eQTL_region.txt \
           >> ${tmp_path}/ubclung/ubclung_lookup_suggestive_july2025.txt
done

head -n1 ${tmp_path}/ubclung/eQTL_region_stat/${CREDSET}_eQTL_region.txt > ${tmp_path}/ubclung/header_ubclung_lookup
cat ${tmp_path}/ubclung/header_ubclung_lookup ${tmp_path}/ubclung/ubclung_lookup_suggestive_july2025.txt \
    > input/ubclung_lookup_suggestive_july2025.txt

#look-up for cs vars:
##manually create the input/pos_credset2025 with only pos info.
grep -F -f input/pos_credset2025 input/ubclung_lookup_suggestive_july2025.txt \
    > /rfs/TobinGroup/nnp5/data/ubclung_lookup_suggestive_july2025.txt


##Run colocalisation:
for c in ${!cs[*]}; do

  N=`cat ${tmp_path}/ubclung/${cs[c]}_UBCLung_probesets.txt | wc -l`

  sbatch --array=1-${N}%20 --export=CREDSET="${cs[c]}" ./src/coloc_UBClung/003_submit_coloc_susie_lung_eQTL.sh

 sleep 5

done

######QUALITY CHECKS:
#to find number of genes: 14
wc -l ${tmp_path}/ubclung/SA_*_UBCLung_probesets.txt | sed 's/_/ /g' | sort -k 3,4 -g | awk '{print $1}'
##Check that all genes for each tissue have been analysed:
grep "UBClung" ${tmp_path}/logerror/coloc_susie_UBCLung*.out | awk -F ":" '{print $1}' | sort -u | wc -l

##Check how many genes per tissue have been analysed for colocalisation:
ls -lthr  ${tmp_path}/results/ubclung/*all_coloc.rds | grep "ubclung" | wc -l
ls -lthr ${tmp_path}/results/ubclung/*all_susie*.rds | grep "UBCLung" | wc -l


#Find statistically significant colocalisation results for UBCLung, and add results into var2gene_raw.xlsx:
#R gave error for xlsx and Java, so I find statistically significant colocalisation results for GTExV8 and eqtlGen eQTL, and add results into var2gene_raw.xlsx:
Rscript ./src/coloc_UBClung/004_concat_coloc_results_ubclung.R

#cp all results into /rfs/:
cp ${tmp_path}/results/coloc*.tsv /rfs/TobinGroup/nnp5/data/


#pQTL look-up:
#check if proteins identified by one study are present in the other ?
awk '{print $13}' /scratch/ukb/kaf19/Noemi_V2G/scallop/scallop_* | awk -F "/|.txt" '{print $6}' | sort -u > ../Thesis_analysis/input/protein_scallop_all
awk '{print $2}' /scratch/ukb/kaf19/Noemi_V2G/ukb_olink/ukb_olink_* | sort -u | grep -v "CHROM" > ../Thesis_analysis/input/protein_ukb_all