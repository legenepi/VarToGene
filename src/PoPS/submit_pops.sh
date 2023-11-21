#!/bin/bash

#SBATCH --job-name=SA_POPS
#SBATCH --output=/scratch/gen1/nnp5/Var_to_Gen_tmp/logerror/%x-%j.out
#SBATCH --time=0:30:0
#SBATCH --mem=50gb
#SBATCH --account=gen1
#SBATCH --export=NONE

#intermediate files in:
tmp_path="/scratch/gen1/nnp5/Var_to_Gen_tmp/"

module load python3

python /data/gen1/LF_HRC_transethnic/PoPS/pops/pops.feature_selection.py \
    --features /data/gen1/LF_HRC_transethnic/PoPS/data/PoPS.features.txt.gz \
    --gene_results ${tmp_path}/pops/SA \
	  --out ${tmp_path}/pops/SA
