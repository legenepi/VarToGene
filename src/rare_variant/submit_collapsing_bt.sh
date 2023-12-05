#!/bin/bash
#Rationale: gene-based collapsing rare variant ExWAS - Based on Nick's script
#guidelines at:
#https://dnanexus.gitbook.io/uk-biobank-rap/science-corner/using-regenie-to-generate-variant-masks
#https://community.dnanexus.com/s/question/0D582000004zTItCAM/burden-test-from-regenie-in-ukb-wes-data-the-result-does-not-contain-alleles-or-beta-and-se
#https://github.com/rgcgithub/regenie/issues/327

dx login
dx select Severe_asthma


dx mkdir /analysis/collapsing_Hg38/
dx cd /analysis/collapsing_Hg38/
dx upload /scratch/gen1/nrgs1/rare_variant/collapsing_Hg38/ukb23158_500k_OQFE.annotations.txt.gz
dx upload /scratch/gen1/nrgs1/rare_variant/collapsing_Hg38/ukb23158_500k_OQFE.sets.txt.gz

#Create the backmask.def file:
printf "M1 LoF\nM3 LoF,missense(5/5)\nM4 LoF,missense(5/5),missense(>=1/5)" > /scratch/gen1/nnp5/Var_to_Gen_tmp/rare_variant/backmask.def
dx upload /scratch/gen1/nnp5/Var_to_Gen_tmp/rare_variant/backmask.def
dx cd ../..

SEED=/analysis/
BASE=/collapsing_Hg38/
EXOME_PATH="/Bulk/Exome sequences/Population level exome OQFE variants, PLINK format - final release"
ANNOT=ukb23158_500k_OQFE.annotations.txt.gz
SETLIST=ukb23158_500k_OQFE.sets.txt.gz
EXCLUDE=ukb23158_500k_OQFE.90pct10dp_qc_variants.txt
MASK=backmask.def
AAF_BINS=0.01,0.001,0.0001,0.00001
JOINT_TESTS=acat
MAXAFF=0.01
TESTS=acatv,skato
THREADS=16

if [ $# -gt 1 ]; then
    BUILD_MASK="--build-mask sum"
    OUT=${OUT}_sum
    NAME=${NAME}_sum
else
    BUILD_MASK="--write-mask --write-mask-snplist"
fi

for CHR in {1..22}
do
  NAME=sa_collapsing_chr${CHR}_gene_p
  BED=ukb23158_c${CHR}_b0_v1
  OUT=${BED}_sa_collapsing_backman_gene_p
  REGENIE_CMD="regenie \
    --bed $BED \
    --exclude $EXCLUDE \
    --extract WES_c${CHR}_snps_qc_pass.snplist \
    --step 2 \
    --pred sa_results_pred.list \
    --phenoFile demo_EUR_pheno_cov_broadasthma_app88144.txt \
    --covarFile demo_EUR_pheno_cov_broadasthma_app88144.txt \
    --phenoCol broad_pheno_1_5_ratio \
    --covarCol age_at_recruitment,age2,PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10,genetic_sex \
    --bt \
    --anno-file ${ANNOT} \
    --set-list ${SETLIST} \
    --mask-def ${MASK} $BUILD_MASK \
    --firth --approx --pThresh 0.01 \
    --bsize 400 \
    --minMAC 3 \
    --aaf-bins $AAF_BINS \
    --joint $JOINT_TESTS \
    --vc-maxAAF $MAXAFF \
    --vc-tests $TESTS \
    --rgc-gene-p \
    --threads=$THREADS \
    --out $OUT"
  dx run swiss-army-knife \
    -iin="${EXOME_PATH}/${BED}.bed" \
    -iin="${EXOME_PATH}/${BED}.bim" \
    -iin="${EXOME_PATH}/${BED}.fam" \
    -iin="${SEED}/WES_c${CHR}_snps_qc_pass.snplist" \
    -iin="/demo_EUR_pheno_cov_broadasthma_app88144.txt" \
    -iin="${SEED}/sa_results_pred.list" \
    -iin="${SEED}/sa_results_1.loco.gz" \
    -iin="${SEED}/${BASE}/${ANNOT}" \
    -iin="${EXOME_PATH}/helper_files/${EXCLUDE}" \
    -iin="${SEED}/${BASE}/${SETLIST}" \
    -iin="${SEED}/${BASE}/${MASK}" \
    -icmd="$REGENIE_CMD" \
    --instance-type "mem1_ssd1_v2_x16" \
    --name="$NAME" \
    --destination="${SEED}/${BASE}" \
    --brief --yes --allow-ssh
done

#Downlaod files in ALICE3:
mkdir /data/gen1/UKBiobank_500K/severe_asthma/Noemi_PhD/data/genecollaps_ExWAS
dx download project-GGzFY70JBJzVx22v4Yj980J1:/analysis/collapsing_Hg38/ukb23158_c*


