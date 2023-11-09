#!/bin/bash

#SBATCH --job-name=lookup_lungeqtl
#SBATCH --output=/scratch/gen1/nnp5/Var_to_Gen_tmp/logerror/%x-%j.out
#SBATCH --time=1:0:0
#SBATCH --mem=30gb
#SBATCH --account=gen1
#SBATCH --export=NONE

module load R

Rscript src/coloc_UBClung/000_run_lookup_lung_eQTL.r "UBCLung" $CREDSET

#To be submitted as:
#for c in ${!cs[*]}; do
#
#  sbatch --export=CREDSET="${cs[c]}" ./src/coloc_UBClung/000_submit_lookup_lung_eQTL.sh
#
#done