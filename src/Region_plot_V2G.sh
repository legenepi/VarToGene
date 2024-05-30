#!/bin/bash

#SBATCH --job-name=regionplot
#SBATCH --output=/scratch/gen1/nnp5/Var_to_Gen_tmp/logerror/%x-%j.out
#SBATCH --time=1:0:0
#SBATCH --mem=30gb
#SBATCH --account=gen1
#SBATCH --export=NONE

#Rationale: visualise V2G results in the different credible sets locus- R2 with regards to SNP with highest average PIP.
module load R

#workdir: Var_to_Gene/
PATH_data="/data/gen1/UKBiobank_500K/severe_asthma/Noemi_PhD/data"
#mkdir /scratch/gen1/nnp5/Var_to_Gen_tmp/regionalplot
PATH_TMP="/scratch/gen1/nnp5/Var_to_Gen_tmp/regionalplot"

##digest FAVOR files:
#awk -F "," '{print $1, $2, $3, $8, $9, $10, $11, $12}' input/FAVOR_credset_chrpos38_2023_08_08.csv \
#    > ${PATH_TMP}/FAVOR_credset_annotations_digest_08_08_23.csv


##Calculate R2 with respect to variant with highest PIP in each locus:
##Use the 10,000 European ancestry individuals in /data/gen1/LF_HRC_transethnic/LD_reference/EUR_UKB:
##find all credible set variants in EUR_UKB by chr and pos (alleles can be flipped):
#awk '{print $3" "$3"_"$4"_"}' /data/gen1/UKBiobank_500K/severe_asthma/Noemi_PhD/data/replsugg_valid_credset.txt | \
#    sed 's/ //g' | \
#    grep -F -f - /data/gen1/LF_HRC_transethnic/LD_reference/EUR_UKB/ukb_imp_chr*_EUR_selected_nodups.bim | \
#    awk '{print $2}' > ${PATH_TMP}/cs_variants_EUR_UKB

##find the SNPs in EUR_UKB by chr and pos (alleles can be flipped):
#awk '{print $2}'input/highest_PIP_sentinels | \
#    awk -F "_" '{print $1, $1"_"$2"_"}' | sed 's/ /   /g' | \
#    grep -F -f - /data/gen1/LF_HRC_transethnic/LD_reference/EUR_UKB/ukb_imp_chr*_EUR_selected_nodups.bim | \
#    awk '{print $2}' > ${PATH_TMP}/highest_PIP_sentinels_EUR_UKB

for line in {1..17}
do
SNP=$(awk -v row="$line" ' NR == row {print $0 } ' ${PATH_TMP}/highest_PIP_sentinels_EUR_UKB)
chr=$(awk -F "_" -v row="$line" ' NR == row {print $1 } ' ${PATH_TMP}/highest_PIP_sentinels_EUR_UKB)
SNP_tmp=$(awk -F "_" -v row="$line" ' NR == row {print $1"_"$2 } ' ${PATH_TMP}/highest_PIP_sentinels_EUR_UKB)
start=$(grep ${SNP_tmp}"_" input/highest_PIP_sentinels | awk '{print $1}' | awk -F "_" '{print $(NF-1)}')
end=$(grep ${SNP_tmp}"_" input/highest_PIP_sentinels | awk '{print $1}' | awk -F "_" '{print $(NF)}')

#create R2 according to leading p-value:
#Calculate R2 with respect to the leading SNP:
#module unload plink2/2.00a
#module load plink
#plink --bfile /data/gen1/LF_HRC_transethnic/LD_reference/EUR_UKB/ukb_imp_chr${chr}_EUR_selected_nodups \
#    --chr ${chr} --from-bp ${start} --to-bp ${end} \
#    --allow-no-sex --r2 inter-chr --ld-snp ${SNP} --ld-window-r2 0 \
#    --out ${PATH_TMP}/${SNP}_${start}_${end}_ld_file

Rscript src/Region_plot_V2G_2.R \
    ${PATH_TMP}/${SNP}_${start}_${end}_ld_file.ld \
    ${chr}_${SNP}_${start}_${end} \
    output/region_plots_V2G/rp_v2g_${chr}_${SNP}_${start}_${end}.pdf \
    ${start} ${end} ${chr} ${SNP}

done