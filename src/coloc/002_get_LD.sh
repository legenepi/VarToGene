#!/bin/bash

#SBATCH --job-name=ld
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4
#SBATCH --output=/scratch/gen1/nnp5/Var_to_Gen_tmp/logerror/%x-%j.out
#SBATCH --time=1:0:0
#SBATCH --mem=24gb
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=nnp5@le.ac.uk
#SBATCH --account=gen1
#SBATCH --export=NONE

module load plink/1.9-beta6.27-vhw5dr2
module load R

#For GTExv8:
pairs_lookup_file="/scratch/gen1/nnp5/Var_to_Gen_tmp/gtex_Pairs_lookup.txt"
DIR="/scratch/gen1/nnp5/Var_to_Gen_tmp"

#For eqtlGen:
#pairs_lookup_file="/scratch/gen1/nnp5/Var_to_Gen_tmp/eqtlgen/eqtlGenWB_Pairs_lookup.txt"
#DIR="/scratch/gen1/nnp5/Var_to_Gen_tmp/eqtlgen"

#For UBCLung:
#pairs_lookup_file="/scratch/gen1/nnp5/Var_to_Gen_tmp/ubclung/UBCLung_Pairs_lookup.txt"
#DIR="/scratch/gen1/nnp5/Var_to_Gen_tmp/ubclung"

cd ${DIR}

while IFS='' read -r LINE || [ -n "${LINE}" ]
do

chr=`echo ${LINE} | awk '{print $4}'`
pheno=`echo ${LINE} | awk '{print $1}'`
signal=`echo ${LINE} | awk '{print $3}'`
eQTL=`echo ${LINE} | awk '{print $6}'`
start=`echo ${LINE} | awk '{print $7}'`
end=`echo ${LINE} | awk '{print $8}'`
anc=`echo ${LINE} | awk '{print $2}'`


name="${DIR}/${signal}_${pheno}_${eQTL}"

plink_script="${name}_plink.sh"
 if [ -f "$plink_script" ]; then
     echo "PLINK script exists, now deleted and created a new one"
       rm $plink_script
 fi

 if [ ${anc}=="EUR" ]; then
    ref_ld="/data/gen1/LF_HRC_transethnic/signal_selection/LDpanel/EUR_UKB/ukb_imp_chr${chr}_EUR_selected"
 else
    ref_ld="/data/gen1/LF_HRC_transethnic/signal_selection/LDpanel/AFR_1000G/ALL.chr${chr}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes_Amin_Amax_noesv_AFR"
 fi

# Create backslash string
bs="\\"

#create regions to extract
region_file="${name}_region.txt"
echo -e ${chr}'\t'${start}'\t'${end} > ${region_file}

# Create script file
cat <<EOT>> ${plink_script}

/home/n/nnp5/software/plink2 $bs
--bfile $ref_ld $bs
--extract range ${region_file} $bs
--make-bed $bs
--out ${name}
EOT

sh ${plink_script}

# Now recode
plink_scriptR="${name}_plink_recode.sh"
if [ -f "$plink_script" ]; then
    echo "PLINK script exists, now deleted and created a new one"
       rm $plink_scriptR
fi

# Create script file
cat <<EOT >> ${plink_scriptR}

/home/n/nnp5/software/plink2 $bs
--bfile ${name} $bs
--out ${name} $bs
--recode A

EOT

sh ${plink_scriptR}

done < ${pairs_lookup_file}