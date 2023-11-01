#!/bin/bash

#SBATCH --job-name=edit_eQTLGen
#SBATCH --output=/scratch/gen1/nnp5/Var_to_Gen_tmp/logerror/%x-%j.out
#SBATCH --time=8:0:0
#SBATCH --mem=80gb
#SBATCH --account=gen1
#SBATCH --export=NONE

module load R

Rscript src/coloc/000_run_edit_eQTLGen.R
