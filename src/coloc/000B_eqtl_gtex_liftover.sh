#!/bin/bash

#==================================================================
# File:        000B_eqtl_gtex_liftover.sh
# Project:     SevereAsthma
# Author:      NNP- edited from KC
# Date:        31 October 2023
# Rationale:   Liftover eQTL coordinates from GRCh37 to GRCh38.
# Submission:  qsub 02B_eqtl_gtex_liftover.sh
#==================================================================

#SBATCH --job-name=eqtl_gtex_liftover
#SBATCH --output=/scratch/gen1/nnp5/Var_to_Gen_tmp/logerror/%x-%j.out
#SBATCH --time=23:0:0
#SBATCH --mem=100gb
#SBATCH --account=gen1
#SBATCH --export=NONE
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=nnp5@le.ac.uk


FILE_PATH="/scratch/gen1/nnp5/Var_to_Gen_tmp/liftover_gtexv8/bed"
liftOver="/home/n/nnp5/software/liftOver"
chain="/data/gen1/reference/liftOver/hg38ToHg19.over.chain"

cd ${FILE_PATH}

for filename in *.hg38.bed
do

  name=$(basename ${filename} .hg38.bed)

  ${liftOver} \
  ${filename} \
  ${chain} \
  ${name}.hg19.bed \
  ${name}.unlifted.bed

done
