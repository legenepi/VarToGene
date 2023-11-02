#!/bin/bash

#SBATCH --job-name=coloc_susie_eqtlgen
#SBATCH --output=/scratch/gen1/nnp5/Var_to_Gen_tmp/logerror/%x-%j.out
#SBATCH --time=0:30:0
#SBATCH --mem=32gb
#SBATCH --account=gen1
#SBATCH --export=NONE

tmp_path="/scratch/gen1/nnp5/Var_to_Gen_tmp/eqtlgen/"

i=$SLURM_ARRAY_TASK_ID
read GENE < <(awk -F'\t' 'NR == '$((i))' { print $1 }' ${tmp_path}${CREDSET}_eqtlGenWB_genes.txt)

module load R

Rscript src/coloc/003_run_coloc_susie_eQTLGen.R "eqtlGenWB" $CREDSET $GENE