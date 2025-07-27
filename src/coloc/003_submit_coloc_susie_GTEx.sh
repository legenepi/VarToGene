#!/bin/bash

#SBATCH --job-name=coloc_susie_gtex
#SBATCH --output=/scratch/gen1/nnp5/Var_to_Gen_tmp/logerror/%x-%j.out
#SBATCH --time=0:30:0
#SBATCH --mem=20gb
#SBATCH --account=gen1
#SBATCH --export=NONE

tmp_path="/scratch/gen1/nnp5/Var_to_Gen_tmp/"


i=$SLURM_ARRAY_TASK_ID
read GENE < <(awk -F'\t' 'NR == '$((i))' { print $1 }' ${tmp_path}${CREDSET}_${TISSUE}_genes.txt)

module load R

Rscript ./src/coloc/003_run_coloc_susie_GTEx.R $TISSUE $CREDSET $GENE





#for c in ${!cs[*]}; do

#  N=`cat ${tmp_path}${cs[c]}_${tissue[t]}_genes.txt | wc -l`

#  sbatch --array=1-${N}%20 --export=TISSUE="${tissue}",CREDSET="${cs[c]}" ${tmp_path}003_submit_coloc_susie_GTEx.sh

# sleep 5

#done