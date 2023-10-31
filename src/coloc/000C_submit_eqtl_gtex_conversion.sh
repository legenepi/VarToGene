#!/bin/bash

#==================================================================
# File:        000A_submit_eqtl_gtex_conversion.sh
# Project:     SevereAsthma
# Author:      NNP
# Date:        31 October 2023
# Rationale:   Submit .
# Submission:  sbatch --array=1-22 src/coloc/000A_submit_eqtl_gtex_conversion.sh
#==================================================================

#SBATCH --job-name=eqtl_gtex_extraction
#SBATCH --output=/scratch/gen1/nnp5/Var_to_Gen_tmp/logerror/%x-%j.out
#SBATCH --time=3:0:0
#SBATCH --mem=100gb
#SBATCH --account=gen1
#SBATCH --export=NONE
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=nnp5@le.ac.uk
#SBATCH --array=1-22           # array ranks to run


i=$SLURM_ARRAY_TASK_ID
module load R
Rscript src/coloc/000C_eqtl_gtex_conversion.R ${i}