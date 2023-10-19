#!/bin/bash

#SBATCH --job-name=SA_eqtl_lookup
#SBATCH --output=/scratch/gen1/nnp5/Var_to_Gen_tmp/logerror/%x-%j.out
#SBATCH --time=0:30:0
#SBATCH --mem=16gb
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=nnp5@le.ac.uk
#SBATCH --account=gen1
#SBATCH --export=NONE

#intermediate files in:
tmp_path="/scratch/gen1/nnp5/Var_to_Gen_tmp/"

module load R
Rscript ./src/coloc/001_run_eqtl_lookup_GTEx.R $TISSUE $CREDSET
