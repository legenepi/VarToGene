#!/bin/bash

#SBATCH --job-name=SA_eqtlGen_lookup
#SBATCH --output=/scratch/gen1/nnp5/Var_to_Gen_tmp/logerror/%x-%j.out
#SBATCH --time=2:00:00
#SBATCH --mem=72gb
#SBATCH --account=gen1
#SBATCH --export=NONE


module load R

Rscript src/coloc/001_run_eqtl_lookup_eQTLGen.R "eqtlGenWB" $CREDSET

#to be submit as:
#for c in ${!cs[@]}; do
#
#    sbatch --export=CREDSET="${cs[c]}" ./src/coloc/001_submit_eqtl_lookup_eQTLGen.sh
#
#done


