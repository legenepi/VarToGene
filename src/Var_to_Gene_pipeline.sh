#!/usr/bin/env bash

#Rationale: pipeline for Variant to Gene mapping analysis. Output from each analysis: list of genes



module load R
#intermediate files in:
tmp_path="/scratch/gen1/nnp5/Var_to_Gen_tmp/"


################
#1.VARIANT ANNOTATION
Rscript src/Variant_annotation_FAVOR.R
################


################
#2.1 COLOCALISATION - eQTL
#GTExV8, eQTLGen, UBCLung
################
##Exclude chromosome 6 - no eQTL colocalisation for this locus.
##Genomic boundaries: PIP-max causal variant +/- 1Mb

cs=('SA_8_81292599_C_A' 'SA_6_90963614_AT_A' 'SA_5_110401872_C_T' 'SA_2_242692858_T_C' 'SA_15_67442596_T_C' 'SA_12_56435504_C_G' 'SA_11_76296671_G_A' 'SA_9_6209697_G_A' 'SA_5_131885240_C_G' 'SA_3_33042712_C_T' 'SA_2_102913642_AAAAC_A' 'SA_17_38168828_G_A' 'SA_16_27359021_C_A'  'SA_15_61068954_C_T' 'SA_12_48202941_T_C' 'SA_10_9064716_C_T')
cs_all="/data/gen1/UKBiobank_500K/severe_asthma/Noemi_PhD/data/replsugg_valid_credset.txt"

###GTExv8 eQTL###

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
#with parameters for GTExV8
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
#wc -l SA_*_${tissue}_genes.txt | sed 's/_/ /g' | sort -k 3,4 -g | awk '{print $1}'
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

#gtex converted in gene symbol:
#https://www.biotools.fr/human/ensembl_symbol_converter
#added in GTExV8_eQTL_genes_symbol table in the input/var2gene.xlsx file.

###eqtlGen eQTL###
mkdir ${tmp_path}/results/eqtlgen
mkdir ${tmp_path}/eqtlgen
dos2unix src/coloc/000_submit_edit_eQTLGen.sh src/coloc/000_run_edit_eQTLGen.R
chmod +x src/coloc/000_submit_edit_eQTLGen.sh src/coloc/000_run_edit_eQTLGen.R
sbatch src/coloc/000_submit_edit_eQTLGen.sh

#Create files for eQTLGen colocalisation:
dos2unix src/coloc/001_submit_eqtl_lookup_eQTLGen.sh src/coloc/001_run_eqtl_lookup_eQTLGen.R
chmod +x src/coloc/001_submit_eqtl_lookup_eQTLGen.sh src/coloc/001_run_eqtl_lookup_eQTLGen.R

for c in ${!cs[@]}; do

    sbatch --export=CREDSET="${cs[c]}" ./src/coloc/001_submit_eqtl_lookup_eQTLGen.sh

done

##Get the LD matrix:
##Create the file with gtex-locus pairs:
dos2unix src/coloc/002_prepare_LDinput_eqtlgen.R
chmod +x src/coloc/002_prepare_LDinput_eqtlgen.R
Rscript src/coloc/002_prepare_LDinput_eqtlgen.R "eqtlGenWB"

##Get LD:
#with parameters for eqtlGen
sbatch src/coloc/002_get_LD.sh
##to see if errors in the job: grep "Error" /scratch/gen1/nnp5/Var_to_Gen_tmp/logerror/ld-170767.out

##Run the colocalisation for GTExV8:
#.sh will run .R script:
mkdir ${tmp_path}/results/eqtlgen
dos2unix src/coloc/003_submit_coloc_susie_eQTLGen.sh
chmod +x src/coloc/003_run_coloc_susie_eQTLGen.R
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


#Find statistically significant colocalisation results for GTExV8 and eqtlGen eQTL, and add results into var2gene_raw.xlsx:
Rscript ./src/coloc/004_concat_coloc_results.R


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
##to see if errors in the job: grep "Error" /scratch/gen1/nnp5/Var_to_Gen_tmp/logerror/ld-295744.out

##Create UBCLung eQTL regional data from eQTL summary stats:
#mkdir ${tmp_path}/ubclung/eQTL_region_stat/
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

  N=`cat ${tmp_path}ubclung/${cs[c]}_UBCLung_probesets.txt | wc -l`

  sbatch --array=1-${N}%20 --export=CREDSET="${cs[c]}" ./src/coloc_UBClung/003_submit_coloc_susie_lung_eQTL.sh

 sleep 5

done

######QUALITY CHECKS:
#to find number of genes: 468
wc -l ${tmp_path}/ubclung/SA_*_UBCLung_probesets.txt | sed 's/_/ /g' | sort -k 3,4 -g | awk '{print $1}'
##Check that all genes for each tissue have been analysed:
grep "UBClung" ${tmp_path}/logerror/coloc_susie_UBCLung*.out | awk -F ":" '{print $1}' | sort -u | wc -l

##Check how many genes per tissue have been analysed for colocalisation:
ls -lthr  ${tmp_path}/results/ubclung/*all_coloc.rds | grep "ubclung" | wc -l
ls -lthr ${tmp_path}/results/ubclung/*all_susie*.rds | grep "UBCLung" | wc -l


#Find statistically significant colocalisation results for UBCLung, and add results into var2gene_raw.xlsx:
#R gave error for xlsx and Java, so I find statistically significant colocalisation results for GTExV8 and eqtlGen eQTL, and add results into var2gene_raw.xlsx:
Rscript ./src/coloc_UBClung/004_concat_coloc_results_ubclung.R

#Added genes into Var_to_Gene/input/var2genes_raw.xlsx file.

#Merge all the gene from eQTL colocalisation:
Rscript src/merge_genes_eqtl_coloc.R

################
#2.2 APPROXIMATE COLOCALISATION - pQTL (LOOK-UP FOR ALL CREDIBLE SET VARIANTS)
#UKB, SCALLOP, deCODE
################
##Exclude chromosome 6 - no eQTL colocalisation for this locus.
##Genomic boundaries: PIP-max causal variant +/- 1Mb

###UKBIOBANK pQTL LOOK-UP###
#Nick generate the pQTL dataset, as for /data/gen1/UKBiobank/olink/pQTL/README.txt:
# OLINK pQTL analysis
#* 48,195 European samples (as defined in Shrine et al. 2023)
#* 1463 proteins (proteins.txt)
#* Phenotype: untransformed log2 fold protein levels (olink_pheno.txt)
#* Covariates: olink batch, age, sex, genotyping array, PC1-10 (olink_covar.txt)
#* 44.8 million MAC >= 5 variants from HRC+UK10K imputation (/data/ukb/imputed_v3)
#* Additive model with regenie (run_step1.sh & run_step2.sh)
#* Full results for each protein tabixed in results directory
#* Script pqtl_lookup.sh for look ups, run with no arguments for usage info

# Nick created a look up script for UKB pQTL:
#/data/gen1/UKBiobank/olink/pQTL/pqtl_lookup.sh -s <RSID> -r <CHR:START-END> [ -f <PROTEINS FILE> ] [ -p <P THRESHOLD> ] [-h]
#UKB pQTL is in GRCh38, need to find sentinel variants position in GRCh38:
mkdir ${tmp_path}/ukb_pqtl
#input for liftOver:
awk -F ' ' 'NR > 1 {print "chr"$3, $4, $4+1}' $cs_all \
    > ${tmp_path}/ukb_pqtl/cs_vars_liftover_input.txt
## download the chain file b37 to b38
wget -P /home/n/nnp5/software/liftover https://hgdownload.soe.ucsc.edu/gbdb/hg19/liftOver/hg19ToHg38.over.chain.gz

/home/n/nnp5/software/liftover/liftOver \
    ${tmp_path}/ukb_pqtl/cs_vars_liftover_input.txt \
    /home/n/nnp5/software/liftover/hg19ToHg38.over.chain.gz \
    ${tmp_path}/ukb_pqtl/cs_vars_liftover_output.bed \
    ${tmp_path}/ukb_pqtl/cs_vars_liftover_unlifted.bed

##divide cs vars into the respective credset with b38 location:
dos2unix src/pQTL_coloc/000_preprocess_cs_b38.R
chmod +x src/pQTL_coloc/000_preprocess_cs_b38.R
Rscript src/pQTL_coloc/000_preprocess_cs_b38.R \
     $cs_all \
     $tmp_path/ukb_pqtl/ \
     ${tmp_path}/ukb_pqtl/cs_vars_liftover_output.bed

#pvalue theshold based on bonferroni correction by the number of measured proteins:
dos2unix src/pQTL_coloc/000_submit_lookup_ukbpqtl.sh
chmod +x src/pQTL_coloc/000_submit_lookup_ukbpqtl.sh
sbatch --array 0-16 src/pQTL_coloc/000_submit_lookup_ukbpqtl.sh

#check that all credible sets variants have been analysed:
find ${tmp_path}/ukb_pqtl/*rs*/ -type f | wc -l

#combine look-up ukb-pqtl files with Nick's script (from /data/gen1/UKBiobank/olink/pQTL/orion_pain):
cd ${tmp_path}/ukb_pqtl/
/home/n/nnp5/PhD/PhD_project/Var_to_Gene/src/pQTL_coloc/001_combine_pqtl.awk *rs*/* \
    > ${tmp_path}/ukb_pqtl/lookup_ukbpqtl.txt
cd /home/n/nnp5/PhD/PhD_project/Var_to_Gene/ #of wherever the folder 'Var_to_Gene' is located
#extract gene showing look-up results:
awk 'NR > 1 {print $2}' ${tmp_path}/ukb_pqtl/lookup_ukbpqtl.txt | sort -u \
    > input/ukbpqtl_var2genes_raw


###deCODE pQTL LOOK-UP###
#Found this script of Nick for decode lookup (from /data/gen1/TSH/coloc_susie/lookup_decode.awk):
#create credible_set.snps: create credible_set.snps in alphabetical order
#($5 < $6 ? $5 : $6)"_"($6 > $5 ? $6 : $5)
#locus85	10_101220474_A_G
#locus85	10_101221275_C_T
mkdir ${tmp_path}/decode_pqtl
#From Chiara/Kayesha script and Nick:
sbatch src/pQTL_coloc/000_submit_lookup_decode.sh
#filter out gene name with significant pQTL:
awk 'NR > 1 && $5 > 0 {print $1}' ${tmp_path}/decode_pqtl/log_pQTL_decode_analysis | sed 's/.txt//g' \
    > input/decode_pqtl_var2genes_raw


###SCALLOP pQTL LOOK-UP###
#From Chiara's script R:\TobinGroup\GWAtraits\Chiara\pQTL_SCALLOP and Nick's script scallop_lookup.awk
mkidr ${tmp_path}/scallop_pqtl
#From Chiara and Nick script:
sbatch src/pQTL_coloc/000_submit_lookup_scallop.sh
#filter out gene name with significant pQTL:
awk 'NR > 1 && $5 > 0 {print $1}' ${tmp_path}/scallop_pqtl/log_pQTL_SCALLOP_analysis | sed 's/.txt//g' \
    > input/scallop_pqtl_var2genes_raw


#Merge genes from the different pQTL look-up analyses:
cat input/ukbpqtl_var2genes_raw input/scallop_pqtl_var2genes_raw input/decode_pqtl_var2genes_raw \
    | awk '{print $2="pQTL", $1}' > input/pqtl_lookup_genes_merged


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
Rscript src/PoPS/PoPS_summary.R

#Added /home/n/nnp5/PhD/PhD_project/Var_to_Gene/input/pops_var2genes_raw.txt to the var2genes_raw.xlsx.

################
#3 NEARBY HUMAN ORTHOLOG MOUSE KO GENE
################
#From Jing's code as well as my own input:
#Download genotype-phenotype data from IMPC (International Mouse Phenotyping consortium)
#https://www.mousephenotype.org/data/release
mkdir ${tmp_path}/mouse_ko
#latest release 2023-07-06:
wget -P ${tmp_path}/mouse_ko/ https://ftp.ebi.ac.uk/pub/databases/impc/all-data-releases/latest/results/genotype-phenotype-assertions-ALL.csv.gz
#latest release 2023-11-22:
wget -P ${tmp_path}/mouse_ko/ http://ftp.ebi.ac.uk/pub/databases/genenames/hcop/human_mouse_hcop_fifteen_column.txt.gz
#run the analysis:
Rscript src/mouse_ko/mouse_ko.r > ${tmp_path}/mouse_ko/output_mouse_ko
#Upload /home/n/nnp5/PhD/PhD_project/Var_to_Gene/input/mouse_ko_genes_raw.txt genes into var2genes_raw.xlsx.

################
#4 NEARBY RARE MENDELIAN DISEASE GENE
################
#From Jing's code as well as my own input:
#Download genotype-phenotype data from https://www.orphadata.com/genes/:
mkdir ${tmp_path}/rare_disease/
#latest release 2023-06-22:
wget -P ${tmp_path}/rare_disease/ https://www.orphadata.com/data/xml/en_product6.xml
#downloaded locally and converted in xlsx in Excel - uploaded in /home/n/nnp5/PhD/PhD_project/Var_to_Gene/input/en_product6.xlsx
#dowload description of the file:
wget -P src/report/ https://www.orphadata.com/docs/OrphadataFreeAccessProductsDescription.pdf
#Downloaded hpo: human phenotype ontology provides a standardized vocabulary of phenotypic abnormalities encounterd in human disease
#latest release June 2023:
#Downloaded locally and converted in xlsx in Excel - uploaded in /home/n/nnp5/PhD/PhD_project/Var_to_Gene/input/en_product4.xlsx
#https://github.com/Orphanet/Orphadata_aggregated/blob/master/Rare%20diseases%20with%20associated%20phenotypes/en_product4.xml
#run the analysis:
Rscript src/rare_disease/rare_disease.r > ${tmp_path}/rare_disease/output_rare_disease
#Upload /home/n/nnp5/PhD/PhD_project/Var_to_Gene/input/rare_disease_genes_raw.txt genes into var2genes_raw.xlsx.


################
#5 RARE VARIANT UKBiobank ANALYSIS in RAP
################
#Prep file ALICE3:
##Match UK Biobank application IDs:
mkdir ${tmp_path}/rare_variant/
Rscript src/rare_variant/000_dataprep_rarevar.R
#Copy the pheno covatiare file in /rfs/:
cp ${tmp_path}/rare_variant/demo_EUR_pheno_cov_broadasthma_app88144.txt /rfs/TobinGroup/data/UKBiobank/application_88144/
#Gene-collapsing rare variant analysis:
#https://github.com/legenepi/rare_collapsing
#Single rare variant analysis: NB - need to do it for rare variant ExWAS ! change filter in plink!
#src/rare_variant/submit_rare_variant.sh