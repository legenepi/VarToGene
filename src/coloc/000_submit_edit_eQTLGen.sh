#!/bin/bash

#SBATCH --job-name=edit_eQTLGen
#SBATCH --output=/scratch/gen1/atw20/pain/log/%x-%j.out
#SBATCH --time=8:0:0
#SBATCH --mem=80gb
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=atw20@le.ac.uk
#SBATCH --account=gen1
#SBATCH --export=NONE

module load R

Rscript /scratch/gen1/atw20/pain/scripts/coloc/000_run_edit_eQTLGen.R
