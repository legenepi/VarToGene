#!/bin/bash

#SBATCH --job-name=ukbpqtol_lookup
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4
#SBATCH --output=/scratch/gen1/nnp5/Var_to_Gen_tmp/logerror/%x-%j.out
#SBATCH --time=72:0:0
#SBATCH --mem=24gb
#SBATCH --account=gen1
#SBATCH --export=NONE


i=$((SLURM_ARRAY_TASK_ID))

#intermediate files in:
tmp_path="/scratch/gen1/nnp5/Var_to_Gen_tmp/"
THRESH=`bc <<< "scale=10; 0.05/1463"`
#locus=('9_rs992969_5709697_6709697' '8_rs7824394_80792599_81792599' '6_rs148639908_90463614_91463614' '5_rs2188962_rs848_131270805_132496500' '5_rs1837253_rs3806932_109901872_110905675' '3_rs35570272_32547662_33547662' '2_rs6761047_242192858_243192858' '2_rs12470864_102426362_103426362' '17_17:38073838_CCG_C_37573838_38573838' '16_rs3024619_26864806_27864806' '15_rs17293632_66942596_67942596' '15_rs11071559_60569988_61569988' '12_rs73107993_47695873_48695873' '12_rs705705_55935504_56935504' '11_rs10160518_75796671_76796671' '10_rs201499805_rs1444789_8542744_9564361' '6_rs9271365_32086794_33086794')
locus=('3_rs778801698_49524027_50524027')
CREDSET="${locus[i]}"
mkdir ${tmp_path}/ukb_pqtl/${CREDSET}

awk 'NR > 1 {print$0}'  ${tmp_path}/ukb_pqtl/SA_${CREDSET}_b38 | while IFS='' read -r LINE || [ -n "${LINE}" ]
do

chr=`echo ${LINE} | awk '{print $3}'`
pos=`echo ${LINE} | awk '{print $11}'`

/data/gen1/UKBiobank/olink/pQTL/pqtl_lookup.sh -r ${chr}:${pos}-${pos} -p ${THRESH} \
    >  ${tmp_path}/ukb_pqtl/${CREDSET}/SA_${chr}_${pos}_ukb_protein_lookup.txt
done

##run job as:
#sbatch --array 0-16 src/pQTL_coloc/000_submit_lookup_ukbpqtl.sh
#for chromsome 3 rs778801698:
#sbatch --array 0 src/pQTL_coloc/000_submit_lookup_ukbpqtl.sh