#!/usr/bin/env bash

#Rationale: pipeline for Variant to Gene mapping analysis for additional credset SNPs - March 2025.
#Output from each analysis: list of genes

module load R
#intermediate files in:
tmp_path="/scratch/gen1/nnp5/Var_to_Gen_tmp"
kath_tmp="/scratch/ukb/kaf19/Noemi_V2G/eqtl_gtex"

################
#1.VARIANT ANNOTATION
################
#From FAVOR website: Two Info:
#Nearest gene for sentinel with highest PIP in each locus
#Functional annotated credset variants according to FANTOM5, ClinVar and Integrative Functional Score criteria
##FAVOR:
##liftover for vars with no rsid:
#liftOver online:
#input/Additional_credset_snps_March2025/Credible_sets.xlsx
#input/Additional_credset_snps_March2025/LiftOver_output_b38_additional_credset_snps.bed
###use FAVOR webtool:
#input/Additional_credset_snps_March2025/
#https://favor.genohub.org/FAVOR_input_additional_credset_SNPs.txt

Rscript src/Variant_annotation_FAVOR_additional_credset_SNPs.R

################
#2.1 COLOCALISATION - eQTL
#KATH RAN LOOK-UP:
#COLOCALISATION TO BE RUN FOR rs705705.12.55935504.56935504 AND rs2188962_rs848.5.131270805.132496500
#rs2188962_rs848.5.131270805.132496500: there are two credset and therefore two PIP sentinel variant - run coloc separate.
#PIP sentinel variant for rs2188962_rs848.5.131270805.132496500: rs1986009 and rs2070729
#CHECK GTExV8 LOOK-UP FOR 5_rs2188962_rs848 BECAUSE IT LOOKS LIKE THE GENES ARE DIFFERENT FROM THE ONES I REPORTED IN
#THE SUPPLEMENTARY TABLE FOR COLOCALISATION RESULTS.
#GTExV8, eQTLGen, UBCLung
################

##variables needed:
tissue=('Stomach' 'Small_Intestine_Terminal_Ileum' 'Lung' 'Esophagus_Muscularis' 'Esophagus_Gastroesophageal_Junction' 'Artery_Tibial' 'Artery_Coronary' 'Artery_Aorta' 'Colon_Transverse' 'Colon_Sigmoid' 'Skin_Sun_Exposed_Lower_leg' 'Skin_Not_Sun_Exposed_Suprapubic')

#Obtain GWASpairs from credible set regions:
#.sh will run .R script:
#credset:
cs=('SA_12_57493727_G_T' 'SA_5_131887986_A_C' 'SA_5_131819921_A_C')
for c in ${!cs[*]}; do
  CREDSET="${cs[c]}"
  Rscript src/coloc/001_run_GWASpairs.R $CREDSET
done

##Get the LD matrix:
##Create the file with gtex-locus pairs:
for t in ${!tissue[*]}; do Rscript ./src/coloc/002_prepare_LDinput.R "${tissue[t]}"; done

##Get LD:
#with parameters for GTExV8
sbatch ./src/coloc/002_get_LD.sh

##to see if errors in the job: grep "Error" /scratch/gen1/nnp5/Var_to_Gen_tmp/logerror/ld-3596057.out

##Run the colocalisation for GTExV8:
#.sh will run .R script:
mkdir ${tmp_path}/results
mkdir ${tmp_path}/results/gtex

#Path for ${cs}_${tissue}_pairs.txt:
#/scratch/ukb/kaf19/Noemi_V2G/eqtl_gtex/tmp/susie_replsugg_credset.rs3024971.12.56993727.57993727_Lung_pairs.txt.gz


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

#create gene file for each tissue:
sort ${kath_tmp}/tmp/susie_replsugg_credset.rs3024971.12.56993727.57993727_${tissue}_genes.txt | uniq \
   > ${tmp_path}/${cs}_${tissue}_genes.txt

sort ${kath_tmp}/tmp/susie_replsugg_credset.rs2188962_rs848.5.131270805.132496500_${tissue}_genes.txt | uniq \
   > ${tmp_path}/${cs}_${tissue}_genes.txt


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

###eqtlGen eQTL###
mkdir ${tmp_path}/results/eqtlgen
mkdir ${tmp_path}/eqtlgen
kath_tmp="/scratch/ukb/kaf19/Noemi_V2G/eQTLgen"
#some data:
#/scratch/ukb/kaf19/Noemi_V2G/eqtl_gtex/tmp/susie_replsugg_credset.rs3024971.12.56993727.57993727_eqtlGenWB_pairs.txt.gz
#/scratch/ukb/kaf19/Noemi_V2G/eqtl_gtex/tmp/susie_replsugg_credset.rs2188962_rs848.5.131270805.132496500_eqtlGenWB_pairs.txt.gz
#/scratch/ukb/kaf19/Noemi_V2G/eqtl_gtex/tmp/susie_replsugg_credset.rs3024971.12.56993727.57993727_eqtlGenWB_genes.txt
#/scratch/ukb/kaf19/Noemi_V2G/eqtl_gtex/tmp/susie_replsugg_credset.rs2188962_rs848.5.131270805.132496500_eqtlGenWB_genes.txt

##Get the LD matrix:
##Create the file with gtex-locus pairs:
Rscript src/coloc/002_prepare_LDinput_eqtlgen.R "eqtlGenWB"

##Get LD:
#with parameters for eqtlGen
sbatch src/coloc/002_get_LD.sh
##to see if errors in the job: grep "Error" /scratch/gen1/nnp5/Var_to_Gen_tmp/logerror/ld-1247283.out

##Run the colocalisation for GTExV8:
#.sh will run .R script:
for c in ${!cs[*]}; do

  N=`cat ${tmp_path}/eqtlgen/${cs[c]}_eqtlGenWB_genes.txt | wc -l`

  sbatch --array=1-${N}%20 --export=CREDSET="${cs[c]}" ./src/coloc/003_submit_coloc_susie_eQTLGen.sh

 sleep 5

done

######QUALITY CHECKS:
#to find number of genes:
wc -l ${tmp_path}/eqtlgen/SA_*_eqtlGenWB_genes.txt | sed 's/_/ /g' | sort -k 3,4 -g | awk '{print $1}'
##Check that all genes for each tissue have been analysed:
grep "eqtlGenWB" ${tmp_path}/logerror/coloc_susie_eqtlgen*.out | awk -F ":" '{print $1}' | sort -u | wc -l

##Check how many genes per tissue have been analysed for colocalisation:
ls -lthr  ${tmp_path}/results/eqtlgen/*all_coloc.rds | grep "eqtlGenWB" | wc -l
ls -lthr ${tmp_path}/results/eqtlgen/*all_susie*.rds | grep "eqtlGenWB" | wc -l

#Find statistically significant colocalisation results for GTExV8 and eqtlGen eQTL, and add results into var2gene_raw.xlsx:
Rscript src/coloc/004_concat_coloc_results_chr3_50024027.R

awk 'NR ==1; $11 == "TRUE" {print $0}' ${tmp_path}/results/coloc_asthma_GTEx.tsv \
    > output/coloc_asthma_GTEx.tsv

awk 'NR ==1; $16 == "TRUE" {print $0}' ${tmp_path}/results/colocsusie_asthma_GTEx.tsv \
    > output/colocsusie_asthma_GTEx.tsv

awk 'NR ==1; $11 == "TRUE" {print $0}' /scratch/gen1/nnp5/Var_to_Gen_tmp/results/coloc_asthma_eqtlgen.tsv \
    > output/coloc_asthma_eqtlgen.tsv

 awk 'NR ==1; $16 == "TRUE" {print $0}' ${tmp_path}/results/colocsusie_asthma_eqtlgen.tsv \
    > output/colocsusie_asthma_eqtlgen.tsv

#gtex converted in gene symbol:
#https://www.biotools.fr/human/ensembl_symbol_converter

#No genes that colocalised with the chose PP.H4.abf threshold.


#Additional for curiosity only:
#If I try filter for PP.H4.abf > 0.85
#
#GTEXv8: There is colocalisation for three genes: ENSG00000004534.14:RBM6, ENSG00000164078.12:MST1R, ENSG00000182179.12:UBA7.
#nsnps PP.H4.abf gene tissue
#1432 0.8969304215521776 ENSG00000004534.14 coronary
#1434 0.8942873733156056 ENSG00000004534.14 not_sun_exposed_suprapubic
#1434 0.8922308534250772 ENSG00000004534.14 sigmoid
#1434 0.8921686387956137 ENSG00000004534.14 muscularis
#1434 0.8921109118112002 ENSG00000004534.14 tibial
#1434 0.8920935039077741 ENSG00000004534.14 transverse
#1434 0.8919127349079411 ENSG00000004534.14 ensg00000004534.14_stomach
#1434 0.8918331099081864 ENSG00000004534.14 sun_exposed_lower_leg
#1434 0.8918261031810106 ENSG00000004534.14 aorta
#1434 0.8916789774689768 ENSG00000004534.14 gastroesophageal_junction
#1434 0.888817774187038 ENSG00000004534.14 ensg00000004534.14_lung
#1434 0.867383436884697 ENSG00000164078.12 muscularis
#1434 0.8562550686662334 ENSG00000164078.12 sun_exposed_lower_leg
#1434 0.8548386547510476 ENSG00000182179.12 aorta
#1434 0.8540520416172054 ENSG00000004534.14 intestine_terminal_ileum
#
#eqtlGen: There is a colocalisation for RBM6:
#nsnps	hit1	hit2	PP.H0.abf	PP.H1.abf	PP.H2.abf	PP.H3.abf	PP.H4.abf	idx1	idx2	snp	n_index	pheno	gene	tissue	coloc_susie
#1384	0	0	3.352464447545198e-4	0.10008140361218355	0.899583349943249	SA	3_50024027_C_CA	RBM6	eqtlGenWB	FALSE

###UBCLung eQTL###
mkdir ${tmp_path}/results/ubclung
mkdir ${tmp_path}/ubclung
dos2unix src/coloc_UBClung/*
chmod +x src/coloc_UBClung/*
for c in ${!cs[*]}; do

  sbatch --export=CREDSET="${cs[c]}" ./src/coloc_UBClung/000_submit_lookup_lung_eQTL.sh

done


##Get the LD matrix:
##Create the file with gtex-locus pairs:
dos2unix src/coloc_UBClung/002_prepare_LDinput_ubclung.R
chmod +x src/coloc_UBClung/002_prepare_LDinput_ubclung.R
Rscript src/coloc_UBClung/002_prepare_LDinput_ubclung.R "UBCLung"

##Get LD:
#with parameters for UBCLung
sbatch src/coloc/002_get_LD.sh
##to see if errors in the job: grep "Error" /scratch/gen1/nnp5/Var_to_Gen_tmp/logerror/ld-1250380.out

##Create UBCLung eQTL regional data from eQTL summary stats:
mkdir ${tmp_path}/ubclung/eQTL_region_stat/
module load tabix
lung_eQTL="/data/gen1/reference/lung_eQTL"

for c in ${!cs[@]}
    do
      CREDSET="${cs[c]}"
      read chr < <(echo $CREDSET | awk -F "_" '{print $2}')
      read pos < <(echo $CREDSET | awk -F "_" '{print $3}')
      pos1=`expr $pos - 1000000`
      pos2=`expr $pos + 1000000`
      if (( pos1 < 0 )); then
          pos1=0
      fi
      tabix -h ${lung_eQTL}/METAANALYSIS_Laval_UBC_Groningen_chr${chr}_formatted.txt.gz ${chr}:${pos1}-${pos2} > ${tmp_path}/ubclung/eQTL_region_stat/${CREDSET}_eQTL_region.txt
done

##Run colocalisation:
dos2unix src/coloc_UBClung/003_submit_coloc_susie_lung_eQTL.sh
dos2unix src/coloc_UBClung/003_run_coloc_susie_lung_eQTL.r
chmod +x src/coloc_UBClung/003_submit_coloc_susie_lung_eQTL.sh
chmod +x src/coloc_UBClung/003_run_coloc_susie_lung_eQTL.r

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

awk 'NR ==1; $13 == "TRUE" {print $0}' ${tmp_path}/results/coloc_asthma_ubclung.tsv \
    > output/coloc_asthma_ubclung_${region}.tsv

awk '$13 == "TRUE" {print $14}' ${tmp_path}/results/coloc_asthma_ubclung.tsv \
    > input/ubclung_gene_${region}.tsv


#Added genes into Var_to_Gene/input/var2genes_raw_chr3_49524027_50524027_rs77880169.xlsx file :
#The only gene for colocalisation among the three method is RBM6 from ubclung.


################
#2.2 APPROXIMATE COLOCALISATION - pQTL (LOOK-UP FOR ALL CREDIBLE SET VARIANTS)
#UKB, SCALLOP, deCODE
################
##Exclude chromosome 6 - no eQTL colocalisation for this locus.
##Genomic boundaries: PIP-max causal variant +/- 1Mb

###UKBIOBANK pQTL LOOK-UP###
#Kath did this

###deCODE pQTL LOOK-UP###
PQTL_PATH="/scratch/gen1/nnp5/Var_to_Gen_tmp/decode_pqtl"
sbatch src/pQTL_coloc/000_submit_lookup_decode.sh
awk 'NR > 1 && $5 > 0 {print $1}' ${PQTL_PATH}/log_pQTL_decode_analysis | sed 's/.txt//g' \
    > /home/n/nnp5/PhD/PhD_project/Var_to_Gene/input/Additional_credset_snps_March2025/decode_pqtl_var2genes_raw_additional_credsetMarch2025
##no results for decode

###SCALLOP pQTL LOOK-UP###
#Kath did this

################
#3 POLYGENIC PRIORITY SCORE (PoPS)
################
#https://github.com/FinucaneLab/pops
#Looking at the github repo and to Jing's code, I compiled PoPS in PoPS.sh and submit_pops.sh:
#PoPS.sh internally calls submit_pops.sh to be submitted as a job - it requires lots of time.
mkdir ${tmp_path}/pops
mkdir ${tmp_path}/pops/results
bash PoPS.sh

#Look at the results and find top score genes within a +/-250Kb from highest-PIP variant for each locus -
#if no top genes in +/- 250Kb, enlarge the window to +/-500Kb:
Rscript src/PoPS/PoPS_summary_additionalcredset_March2025.R
#Put genes in the Var_to_Gene/input/Additional_credset_snps_March2025/var2genes_raw_additionalcredset_March2025.xlsx
cp /scratch/gen1/nnp5/Var_to_Gen_tmp/pops/results/additional_credsetMarch2025_all_results_merged_table.txt \
    input/Additional_credset_snps_March2025/pops_var2genes_additional_credsetMarch2025_all_results_merged_table.txt

################
#3 NEARBY HUMAN ORTHOLOG MOUSE KO GENE
################
#From Jing's code as well as my own input:
#Download genotype-phenotype data from IMPC (International Mouse Phenotyping consortium)
#https://www.mousephenotype.org/data/release
mkdir ${tmp_path}/mouse_ko
#release 2024-12-18 for this additional credset variants:
wget -P ${tmp_path}/mouse_ko/ https://ftp.ebi.ac.uk/pub/databases/impc/all-data-releases/latest/results/genotype-phenotype-assertions-ALL.csv.gz
#release 2024-07-26 for this additional credset variants:
wget -P ${tmp_path}/mouse_ko/ http://ftp.ebi.ac.uk/pub/databases/genenames/hcop/human_mouse_hcop_fifteen_column.txt.gz
#run the analysis:
Rscript src/mouse_ko/mouse_ko_additionalcredset_March2025.r
#Upload the genes into input/Additional_credset_snps_March2025/var2genes_raw_additionalcredset_March2025.xlsx

################
#4 NEARBY RARE MENDELIAN DISEASE GENE
################
#From Jing's code as well as my own input:
#Download genotype-phenotype data from https://www.orphadata.com/genes/:
mkdir ${tmp_path}/rare_disease/
#use the release of 2023-06-22:
#downloaded locally and converted in xlsx in Excel - uploaded in /home/n/nnp5/PhD/PhD_project/Var_to_Gene/input/en_product6.xlsx
#dowload description of the file:
wget -P src/report/ https://www.orphadata.com/docs/OrphadataFreeAccessProductsDescription.pdf
#Downloaded hpo: human phenotype ontology provides a standardized vocabulary of phenotypic abnormalities encounterd in human disease
#use the release of June 2023:
#Downloaded locally and converted in xlsx in Excel - uploaded in /home/n/nnp5/PhD/PhD_project/Var_to_Gene/input/en_product4.xlsx
#https://github.com/Orphanet/Orphadata_aggregated/blob/master/Rare%20diseases%20with%20associated%20phenotypes/en_product4.xml
#run the analysis:
Rscript src/rare_disease/rare_disease_additionalcredset_March2025.r
#Upload genes input/Additional_credset_snps_March2025/var2genes_raw_additionalcredset_March2025.xlsx

################
#5 RARE VARIANT UKBiobank ANALYSIS in RAP
################
#Significant results not in these region !

################
#5 Combine all results for locus and gene - supplementary table with SNPs of interest.
################