#!/usr/bin/env bash

#Rationale: pipeline for Variant to Gene mapping analysis for chromosome 3, region 49524027_50524027, rs778801698.
#Output from each analysis: list of genes

module load R
#intermediate files in:
tmp_path="/scratch/gen1/nnp5/Var_to_Gen_tmp"

region="chr3_49524027_50524027_rs778801698"

################
#1.VARIANT ANNOTATION
################
#From FAVOR website: Two Info:
#Nearest gene for sentinel with highest PIP in each locus
#Functional annotated credset variants according to FANTOM5, ClinVar and Integrative Functional Score criteria
##FAVOR:
##liftover for vars with no rsid:
awk -F "\t" '$1 == "3_rs778801698_49524027_50524027" {print "chr"$3,$4,$4+1}' /data/gen1/UKBiobank_500K/severe_asthma/Noemi_PhD/data/replsugg_valid_credset_chr3_noMHC.txt \
    > input/replsugg_valid_credset_input_liftover_b37_chr3

#liftOver online:
#input/hglft_genome_credset_vars_chr3_49524027_50524027_rs778801698.bed


awk -F "\t" '$1 == "3_rs778801698_49524027_50524027" {print $5"-"$6}' \
    /data/gen1/UKBiobank_500K/severe_asthma/Noemi_PhD/data/replsugg_valid_credset_chr3_noMHC.txt \
    > /scratch/gen1/nnp5/Fine_mapping/tmp_data/alleles_for_favor_input_${region}

awk '{print $1"-"$2}' input/hglft_genome_credset_vars_${region}.bed | sed 's/chr//g' | \
    paste -d "-" - /scratch/gen1/nnp5/Fine_mapping/tmp_data/alleles_for_favor_input_${region} \
    > input/replsugg_credset_chrpos38_${region}.txt

###use FAVOR webtool:
#https://favor.genohub.org/

Rscript src/Variant_annotation_FAVOR_chr3.R \
    input/FAVOR_credset_chrpos38_2024_05_14_${region}.txt.csv ${region}

#copy and paste the gene for FANTOM5-ClinVar-Integrative Functional Score in the varannot_gene sheet of input/var2genes_raw_chr3_49524027_50524027_rs778801698.xlsx"

################
#2.1 COLOCALISATION - eQTL
#GTExV8, eQTLGen, UBCLung
################
##Exclude chromosome 6 - no eQTL colocalisation for this locus.
##Genomic boundaries: PIP-max causal variant +/- 1Mb

cs=('SA_3_50024027_C_CA')
cs_all="/data/gen1/UKBi/replsugg_valid_credset_chr3_noMHC.txt"

###GTExv8 eQTL###
##Create eQTl files in hg19 for Colon Transverse, Colon Sigmoid, Skin_Not_Sun_Exposed_Suprapubic, Skin_Sun_Exposed_Lower_leg
#From GTExV8 .parquet files and in hg38.
mkdir ${tmp_path}/liftover_gtexv8
mkdir ${tmp_path}/liftover_gtexv8/bed

dos2unix src/coloc/000A_eqtl_gtex_extraction.R src/coloc/000B_eqtl_gtex_liftover.sh \
    src/coloc/000C_eqtl_gtex_conversion.R src/coloc/000A_submit_eqtl_gtex_extraction.sh src/coloc/000C_submit_eqtl_gtex_conversion.sh
chmod o+x src/coloc/000A_eqtl_gtex_extraction.R src/coloc/000B_eqtl_gtex_liftover.sh \
    src/coloc/000C_eqtl_gtex_conversion.R src/coloc/000A_submit_eqtl_gtex_extraction.sh src/coloc/000C_submit_eqtl_gtex_conversion.sh

sbatch --array=3 src/coloc/000A_submit_eqtl_gtex_extraction.sh

sbatch src/coloc/000B_eqtl_gtex_liftover.sh

sbatch --array=3 src/coloc/000C_submit_eqtl_gtex_conversion.sh

##variables needed:
tissue=('Stomach' 'Small_Intestine_Terminal_Ileum' 'Lung' 'Esophagus_Muscularis' 'Esophagus_Gastroesophageal_Junction' 'Artery_Tibial' 'Artery_Coronary' 'Artery_Aorta' 'Colon_Transverse' 'Colon_Sigmoid' 'Skin_Sun_Exposed_Lower_leg' 'Skin_Not_Sun_Exposed_Suprapubic')

##Divide the credible sets into separate files:
Rscript ./src/coloc/000_preprocess_cs.R $cs_all $tmp_path/

#Obtain GWASpairs from credible set regions:
#.sh will run .R script:
for c in ${!cs[*]}; do

  sbatch --export=CREDSET="${cs[c]}" ./src/coloc/001_submit_GWASpairs.sh

done

#Create files for GTExV8:
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

##to see if errors in the job: grep "Error" /scratch/gen1/nnp5/Var_to_Gen_tmp/logerror/ld-1236385.out

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
ls -lthr  ${tmp_path}/results/gtex/*all_coloc.rds | grep ${tissue} | wc -l
ls -lthr ${tmp_path}/results/gtex/*all_susie*.rds | grep ${tissue} | wc -l
done

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
##to see if errors in the job: grep "Error" /scratch/gen1/nnp5/Var_to_Gen_tmp/logerror/ld-1247283.out

##Run the colocalisation for GTExV8:
#.sh will run .R script:
mkdir ${tmp_path}/results/eqtlgen
dos2unix src/coloc/003_submit_coloc_susie_eQTLGen.sh
chmod +x src/coloc/003_run_coloc_susie_eQTLGen.R
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