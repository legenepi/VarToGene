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
cs_all="/data/gen1/UKBiobank_500K/severe_asthma/Noemi_PhD/data/replsugg_valid_credset_chr3_noMHC.txt"

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
sbatch --array 0 src/pQTL_coloc/000_submit_lookup_ukbpqtl.sh

#combine look-up ukb-pqtl files with Nick's script (from /data/gen1/UKBiobank/olink/pQTL/orion_pain):
cd ${tmp_path}/ukb_pqtl/
/home/n/nnp5/PhD/PhD_project/Var_to_Gene/src/pQTL_coloc/001_combine_pqtl.awk *rs*/* \
    > ${tmp_path}/ukb_pqtl/lookup_ukbpqtl_${region}.txt
cd /home/n/nnp5/PhD/PhD_project/Var_to_Gene/
cp ${tmp_path}/ukb_pqtl/lookup_ukbpqtl_${region}.txt output/lookup_ukbpqtl_${region}.txt
 #of wherever the folder 'Var_to_Gene' is located
#extract gene showing look-up results:
awk 'NR > 1 {print $2}' ${tmp_path}/ukb_pqtl/lookup_ukbpqtl_${region}.txt | sort -u \
    > input/ukbpqtl_var2genes_raw_${region}
#5 genes

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
    > input/decode_pqtl_var2genes_raw_${region}
#no results

###SCALLOP pQTL LOOK-UP###
#From Chiara's script R:\TobinGroup\GWAtraits\Chiara\pQTL_SCALLOP and Nick's script scallop_lookup.awk
mkidr ${tmp_path}/scallop_pqtl
#From Chiara and Nick script:
sbatch src/pQTL_coloc/000_submit_lookup_scallop.sh
#filter out gene name with significant pQTL:
awk 'NR > 1 && $5 == 1 {print $1}' ${tmp_path}/scallop_pqtl/log_pQTL_SCALLOP_analysis | sed 's/.txt//g' \
    > input/scallop_pqtl_var2genes_raw_${region}
#no results


#Merge genes from the different pQTL look-up analyses:
cat input/ukbpqtl_var2genes_raw_${region} input/scallop_pqtl_var2genes_raw_${region} input/decode_pqtl_var2genes_raw_${region} \
    | awk '{print $2="pQTL", $1}' > input/pqtl_lookup_genes_merged_${region}


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


################
#3 NEARBY HUMAN ORTHOLOG MOUSE KO GENE
################
#From Jing's code as well as my own input:
#Download genotype-phenotype data from IMPC (International Mouse Phenotyping consortium)
#https://www.mousephenotype.org/data/release
mkdir ${tmp_path}/mouse_ko
#release 2023-07-06; release 2024-05-07 for chromosome 3:
wget -P ${tmp_path}/mouse_ko/ https://ftp.ebi.ac.uk/pub/databases/impc/all-data-releases/latest/results/genotype-phenotype-assertions-ALL.csv.gz
#release 2023-11-22; release 2024-05-15 for chromosome 3:
wget -P ${tmp_path}/mouse_ko/ http://ftp.ebi.ac.uk/pub/databases/genenames/hcop/human_mouse_hcop_fifteen_column.txt.gz
#run the analysis:
Rscript src/mouse_ko/mouse_ko.r > ${tmp_path}/mouse_ko/output_mouse_ko_${region}
#Upload /home/n/nnp5/PhD/PhD_project/Var_to_Gene/input/mouse_ko_genes_raw_chr3_49524027_50524027_rs778801698.txt genes into var2genes_raw${region}.xlsx.
cp ${tmp_path}/mouse_ko/results_chr3_noMHC_500kb.csv input/mko_results_500kb${region}.csv
grep "3_rs778801698_49524027_50524027" input/mko_results_500kb${region}.csv | awk -F ","  '{print $6}' | sort -u > input/mko_results_500kb_only_${region}_genes.csv

################
#4 NEARBY RARE MENDELIAN DISEASE GENE
################
#From Jing's code as well as my own input:
#Download genotype-phenotype data from https://www.orphadata.com/genes/:
mkdir ${tmp_path}/rare_disease/
#latest release 2023-06-22:
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
grep "3_" ${tmp_path}/rare_disease/results_chr3_noMHC_500kb_by_gene.csv > input/raredis_results_500kb_only_${region}_genes.csv

################
#5 RARE VARIANT UKBiobank ANALYSIS in RAP
################
#Not on chromosome 3 and in that region !

################
#6 TABLES FOR EACH ANALYSIS AND MERGE GENES FOR GENE PRIORITISATION AND VISUALISATION - add genes from chr 3 rs778801698
################
#TO NOTE: ONLY FOR
#nearest gene
#functional annotation
#eQTL
#mouse_KO
#PoPS
#pQTL
#rare_disease
#SINGLE AND GENE-BASED COLLAPSING ANALYSIS: no genes for variant-to-gene mapping, so I did not add them to this table
Rscript src/Locus_to_genes_table_chr3_rs778801698.R
Rscript src/genes_heatmap.R
cp src/report/var2gene_full_3_rs778801698_49524027_50524027.xlsx /data/gen1/UKBiobank_500K/severe_asthma/Noemi_PhD/data/
cp output/V2G_heatmap_subplots_chr3_noMHC.png /data/gen1/UKBiobank_500K/severe_asthma/Noemi_PhD/data/


################
#7 REGION PLOTS WITH VAR2GENE RESULTS
################
#bash Region_plot_V2G.sh
#Region_plot_V2G.R
#Region_plot_V2G_2.R