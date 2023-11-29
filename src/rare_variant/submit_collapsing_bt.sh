#!/bin/bash

dx mkdir /analysis/collapsing_Hg38/
dx cd /analysis/collapsing_Hg38/
dx upload /scratch/gen1/nrgs1/rare_variant/collapsing_Hg38/ukb23158_500k_OQFE.annotations.txt.gz
dx upload /scratch/gen1/nrgs1/rare_variant/collapsing_Hg38/ukb23158_500k_OQFE.sets.txt.gz

#Create the backmask.def file:
printf "M1 LoF\nM3 LoF,missense(5/5)\nM4 LoF,missense(5/5),missense(>=1/5)" > /scratch/gen1/nnp5/Var_to_Gen_tmp/rare_variant/backmask.def
dx upload /scratch/gen1/nnp5/Var_to_Gen_tmp/rare_variant/backmask.def
dx cd ../..

CHR=$1
BASE=/rare_variants/collapsing_Hg38
REGENIE=regenie_v3.3.gz_x86_64_Centos7_mkl
BED=ukb23158_c${CHR}_b0_v1
EXOME_PATH="/Bulk/Exome sequences/Population level exome OQFE variants, PLINK format - final release"
COVAR="/demo_EUR_pheno_cov_broadasthma_app88144.txt"
PHENO="/demo_EUR_pheno_cov_broadasthma_app88144.txt"
PRED=sa_fit_bt_out_Hg38
ANNOT=ukb23158_500k_OQFE.annotations.txt.gz
SETLIST=ukb23158_500k_OQFE.sets.txt.gz
EXCLUDE=ukb23158_500k_OQFE.90pct10dp_qc_variants.txt
MASK=backmask.def
AAF_BINS=0.01,0.001,0.0001,0.00001
JOINT_TESTS=acat
MAXAFF=0.01
TESTS=acatv,skato
OUT=${BED}_sa_collapsing_backman_gene_p
NAME=sa_collapsing_chr${CHR}_gene_p
THREADS=16

if [ $# -gt 1 ]; then
    BUILD_MASK="--build-mask sum"
    OUT=${OUT}_sum
    NAME=${NAME}_sum
else
    BUILD_MASK="--write-mask --write-mask-snplist"
fi

REGENIE_CMD="cp ${REGENIE} /tmp; chmod u+x /tmp/${REGENIE}; /tmp/${REGENIE} --bed $BED --exclude $EXCLUDE --step 2 --pred ${PRED}_pred.list --covarFile $COVAR --phenoFile $PHENO --bt --anno-file ${ANNOT} --set-list ${SETLIST} --mask-def ${MASK} $BUILD_MASK --firth --approx --pThresh 0.01 --bsize 400 --minMAC 3 --aaf-bins $AAF_BINS --joint $JOINT_TESTS --vc-maxAAF $MAXAFF --vc-tests $TESTS --rgc-gene-p --threads=$THREADS --out $OUT"

dx run swiss-army-knife \
    -iin="${EXOME_PATH}/${BED}.bed" \
    -iin="`[ $CHR != "Y" ] && echo ${EXOME_PATH}/${BED}.bim || echo /rare_variants/${BED}.bim`" \
    -iin="${EXOME_PATH}/${BED}.fam" \
    -iin="${BASE}/$COVAR" \
    -iin="${BASE}/$PHENO" \
    -iin="${BASE}/${PRED}_pred.list" \
    -iin="${BASE}/${PRED}_1.loco" \
    -iin="${BASE}/${PRED}_2.loco" \
    -iin="${BASE}/${PRED}_3.loco" \
    -iin="${BASE}/${PRED}_4.loco" \
    -iin="${BASE}/${PRED}_5.loco" \
    -iin="${BASE}/${PRED}_6.loco" \
    -iin="${BASE}/${PRED}_7.loco" \
    -iin="${BASE}/${PRED}_8.loco" \
    -iin="${BASE}/${PRED}_9.loco" \
    -iin="${BASE}/${PRED}_10.loco" \
    -iin="${BASE}/${PRED}_11.loco" \
    -iin="${BASE}/${PRED}_12.loco" \
    -iin="${BASE}/${PRED}_13.loco" \
    -iin="${BASE}/${PRED}_14.loco" \
    -iin="${BASE}/${PRED}_15.loco" \
    -iin="${BASE}/${PRED}_16.loco" \
    -iin="${BASE}/${PRED}_17.loco" \
    -iin="${BASE}/${PRED}_18.loco" \
    -iin="${BASE}/${PRED}_19.loco" \
    -iin="${BASE}/${PRED}_20.loco" \
    -iin="${BASE}/${PRED}_21.loco" \
    -iin="${BASE}/${PRED}_22.loco" \
    -iin="${BASE}/${PRED}_23.loco" \
    -iin="${BASE}/${ANNOT}" \
    -iin="${EXOME_PATH}/helper_files/${EXCLUDE}" \
    -iin="${BASE}/${SETLIST}" \
    -iin="${BASE}/${MASK}" \
    -iin="${BASE}/${REGENIE}" \
    -icmd="$REGENIE_CMD" \
    --instance-type "mem1_ssd1_v2_x16" \
    --name="$NAME" \
    --destination="$BASE" \
    --brief --yes --allow-ssh
#    -imount_inputs="true" \
