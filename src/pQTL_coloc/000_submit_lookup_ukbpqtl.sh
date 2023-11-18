#!/bin/bash

#SBATCH --job-name=ukbpqtol_lookup
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4
#SBATCH --output=/scratch/gen1/nnp5/Var_to_Gen_tmp/logerror/%x-%j.out
#SBATCH --time=1:0:0
#SBATCH --mem=24gb
#SBATCH --account=gen1
#SBATCH --export=NONE

#intermediate files in:
tmp_path="/scratch/gen1/nnp5/Var_to_Gen_tmp/"
THRESH=`bc <<< "scale=10; 0.05/1463"`

while IFS='' read -r LINE || [ -n "${LINE}" ]
do

chr=`echo ${LINE} | awk -F 'chr| ' '{print $2}'`
pos=`echo ${LINE} | awk '{print $3}'`
start=`echo ${LINE} | awk '{print $3-1000000}'`
end=`echo ${LINE} | awk '{print $3+1000000}'`

echo $chr:$start-$end

/data/gen1/UKBiobank/olink/pQTL/pqtl_lookup.sh -r ${chr}:${start}-${end} -p ${THRESH} \
    >  ${tmp_path}/ukb_pqtl/SA_${chr}_${pos}_ukb_protein_lookup.txt
done < ${tmp_path}/ukb_pqtl/cs_sentinel_vars_liftover_input.txt