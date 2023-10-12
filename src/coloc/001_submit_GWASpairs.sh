#!/bin/bash

#SBATCH --job-name=SA_GWASpairs
#SBATCH --output=/scratch/gen1/nnp5/Var_to_Gen_tmp/%x-%j.out
#SBATCH --time=0:30:0
#SBATCH --mem=16gb
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=nnp5@le.ac.uk
#SBATCH --account=gen1
#SBATCH --export=NONE

module load R

Rscript src/001_run_GWASpairs.R $CREDSET


# Submit jobs

cs=('SA_8_81292599_A_C' 'SA_6_90963614_A_AT' 'SA_5_110401872_T_C' 'SA_2_242692858_C_T' 'SA_15_67442596_C_T' 'SA_12_56435504_G_C' 'SA_11_76296671_A_G' 'SA_9_6209697_A_G' 'SA_6_32586794_T_G' 'SA_5_131885240_G_C' 'SA_3_33042712_T_C' 'SA_2_102913642_A_AAAAC' 'SA_17_38168828_A_G' 'SA_16_27359021_A_C'  'SA_15_61068954_T_C' 'SA_12_48202941_C_T' 'SA_10_9064716_T_C')

module load R

for c in ${!cs[@]}
do
  CREDSET="${cs[c]}"
  Rscript ./src/001_run_GWASpairs.R $CREDSET
done