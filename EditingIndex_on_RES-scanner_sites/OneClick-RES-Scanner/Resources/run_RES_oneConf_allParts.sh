#!/usr/bin/env bash

#for ruunnig RES's all parts for one config


ARGPARSE_DESCRIPTION="ruunnig RES's all parts for one config"      # this is optional
source ~/argparse-bash/argparse.bash || exit 1
argparse "$@" <<EOF || exit 1
parser.add_argument('-o', '--out_dir', type=str, help='Full path of existing output directory', required=True)
parser.add_argument('-c', '--conf_file', type=str, help='Path of config_file', required=True)
parser.add_argument('-r', '--ref_file', type=str, help='Path of genome refernce', required=True)
parser.add_argument('-a', '--additional_arguments', type=str, help='additional argumets and parameters default --ss 0 --editLevel 0.01', required=False, default='--ss 0 --editLevel 0.01')

parser.add_argument('-b', '--bwa_path',type=str, help='Path of BWA ', required=False, default='/private/common/Software/BWA/bwa-0.6.2/bwa')
parser.add_argument('-s', '--samtools_path',type=str, help='Path of samtools ', required=False, default='/private/apps/bin/samtools')
parser.add_argument('-l', '--blat_path',type=str, help='Path of blat ', required=False, default='/private/apps/bin/blat')
EOF

cd $OUT_DIR

SCRIPT_DIR="$(dirname "$(readlink -f "$0")")"
res_S_DIR=${SCRIPT_DIR}/../RES-Scanner/RES-Scanner

perl ${res_S_DIR}/RES-Scanner_alignment/RES-Scanner_alignment.part1.pl --ref $REF_FILE --outDir $OUT_DIR --bwa $BWA_PATH --config $CONF_FILE 

sh "/home/alu/twerski/Scripts/BASH/RES_STEPS/run_step_0_1.sh" $OUT_DIR 

perl ${res_S_DIR}/RES-Scanner_alignment/RES-Scanner_alignment.part2.pl --config $CONF_FILE --ref $REF_FILE --outDir $OUT_DIR -t 12 --bwa $BWA_PATH --samtools $SAMTOOLS_PATH --index 1 

${SCRIPT_DIR}/RES_STEPS/run_step_2_3.sh $OUT_DIR

IDENTIFICATION_OUT_DIR=${OUT_DIR}/"identification/"

mkdir -p $IDENTIFICATION_OUT_DIR;

cd $IDENTIFICATION_OUT_DIR

perl ${res_S_DIR}/RES-Scanner_identification/RES-Scanner_identification.pl --config ${OUT_DIR}"RES_Scanner_indentification_config.txt" --out $IDENTIFICATION_OUT_DIR --genome $REF_FILE --blat $BLAT_PATH --samtools $SAMTOOLS_PATH $ADDITIONAL_ARGUMENTS

sh ${SCRIPT_DIR}/RES_STEPS/run_identif_step1_2_3_4.sh $IDENTIFICATION_OUT_DIR