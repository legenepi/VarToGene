#!/bin/bash

#SBATCH --job-name=pain_eqtl_lookup
#SBATCH --output=/scratch/gen1/atw20/pain/log/%x-%j.out
#SBATCH --time=2:0:0
#SBATCH --mem=72gb
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=atw20@le.ac.uk
#SBATCH --account=gen1
#SBATCH --export=NONE

module load R

Rscript /scratch/gen1/atw20/pain/scripts/coloc/001_run_eqtl_lookup_eQTLGen.R $CREDSET
