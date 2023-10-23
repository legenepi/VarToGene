#!/bin/bash

#SBATCH --job-name=liftover_gtex
#SBATCH --output=/scratch/gen1/nnp5/Var_to_Gen_tmp/logerror/%x-%j.out
#SBATCH --time=0:30:0
#SBATCH --mem=32gb
#SBATCH --account=gen1
#SBATCH --export=NONE

#Script rights: from Dr.Chiara Batini
#Rationale: convert in b37 (/hg19) the eQTL data.
#This script will produce data for all the eQTL tissue. I am interested in Colon_Sigmoid and Colon_Transverse.

tmp_path="/scratch/gen1/nnp5/Var_to_Gen_tmp/"


#for p in $(seq 1 2 49); do sbatch 000_liftover_b38_to_b37_GTExv8.sh  --export=start=${p},end=$((${p}+1)); done

#define liftover tool and chain
liftover="/scratch/gen1/cb334/liftOver"
chain="/data/gen1/reference/liftOver/hg38ToHg19.over.chain"

#define original dataset and output directory
gtex_v8="/data/gen1/reference/GTEx_Analysis_v8_QTLs/GTEx_Analysis_v8_eQTL_all_associations"
#gtex_v8_b37="/data/gen1/reference/GTEx_Analysis_v8_QTLs/GTEx_Analysis_v8_eQTL_all_associations_b37"
output_dir="/scratch/gen1/cb334/liftover_gtexv8"

cd $output_dir

for file in $(ls ${gtex_v8}/*.allpairs.txt.gz | sed -n ${start},${end}p)
do
filename=$(basename ${file} | sed 's/\.allpairs.txt.gz//')

#create b38 bed file
zcat ${gtex_v8}/${filename}.allpairs.txt.gz |
awk 'NR>1' | \
sed 's/_/\t/g' | \
awk {'print $2"        "$3"        "$3+1"        "$2"_"$3"_"$4"_"$5"_"$6'} > ${output_dir}/${filename}_b38.bed

#perform b38 to b37 liftover
$liftover ${output_dir}/${filename}_b38.bed $chain ${output_dir}/${filename}_b37.bed ${output_dir}/${filename}_b38_unlifted.bed


#create a map file b37 b38
echo -e "variant_id_b37\tvariant_id_b38" > ${output_dir}/${filename}_b37_b38_varids

sed 's/_/\t/g' ${output_dir}/${filename}_b37.bed | \
awk '$1==$4 {print $0}' | \
awk ' { print $1"_"$2"_"$6"_"$7"_b37\t"$4"_"$5"_"$6"_"$7"_b38"} ' >> ${output_dir}/${filename}_b37_b38_varids

#recreate qtl files with b37 variant id
awk -v OFS='\t' '
#for the first input file, the map, create the key-value relation (key=b38, value=b37)
NR == FNR {
    map[$2] = $1
    next
}

#for the second input, print the header line
FNR == 1 {
    print
}

# and then for the variants in the first input, substitute the key with the value
FNR > 1 && $2 in map{


	$2 = map[$2]
	print

}' ${output_dir}/${filename}_b37_b38_varids <(zcat ${gtex_v8}/${filename}.allpairs.txt.gz) > ${output_dir}/${filename}.allpairs.b37.txt


gzip -f ${output_dir}/${filename}.allpairs.b37.txt

done