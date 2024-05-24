#!/bin/bash

#SBATCH --job-name=SA_POPS
#SBATCH --output=/scratch/gen1/nnp5/Var_to_Gen_tmp/logerror/%x-%j.out
#SBATCH --time=60:00:00
#SBATCH --mem=50gb
#SBATCH --account=gen1
#SBATCH --export=NONE

#Rationale: PoPS feature selection and score calculation

tmp_path="/scratch/gen1/nnp5/Var_to_Gen_tmp"
popsdir="/data/gen1/LF_HRC_transethnic/PoPS"

module load python3

python ${popsdir}/pops/pops.feature_selection.py \
    --features ${popsdir}/data/PoPS.features.txt.gz \
    --gene_results ${tmp_path}/pops/SA \
	  --out ${tmp_path}/pops/SA \

#for i in {1..22}
#for chr3 rs778801698
i=3
do
    python ${popsdir}/pops/pops.predict_scores.py \
	    --gene_loc ${popsdir}/data/gene_loc.txt \
	    --gene_results ${tmp_path}/pops/SA \
	    --features ${popsdir}/data/PoPS.features.txt.gz \
	    --selected_features  ${tmp_path}/pops/SA.features \
	    --control_features ${popsdir}/data/control.features \
	    --chromosome ${i} \
	    --out ${tmp_path}/pops/SA
done
