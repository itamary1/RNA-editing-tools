#!/usr/bin/env bash


ARGPARSE_DESCRIPTION="ruunnig  RES scanner and then editing index on its results"      
source ~/argparse-bash/argparse.bash || exit 1
argparse "$@" <<EOF || exit 1
parser.add_argument('-m', '--max_workers', type=str, help='number of maximum samples to run in parallel', required=False, default='5')
parser.add_argument('-d', '--project_dir', type=str, help='Path of project directory', required=True)
parser.add_argument('-f', '--fastq_dir', type=str, help='Path of fastq dir', required=True)
parser.add_argument('-su', '--fastq_suffix',type=str, help='fastq files suffix', required=False, default='fastq')
parser.add_argument('-dn', '--dna_file', type=str, help='full Path of dna file', required=True)
parser.add_argument('-ln', '--dna_length',type=str, help='length of reads of DNA', required=False, default='82')
parser.add_argument('-lr', '--rna_length',type=str, help='length of reads of RNA', required=False, default='81')
parser.add_argument('-r', '--ref_file', type=str, help='Path of genome refernce', required=True)
parser.add_argument('-a', '--additional_arguments', type=str, help='additional argumets and parameters default --ss 0 --editLevel 0.01', required=False, default='--ss 0 --editLevel 0.01')
parser.add_argument('-b', '--bwa_dir',type=str, help='Path of BWA ', required=False, default='bwa')
parser.add_argument('-s', '--samtools_dir',type=str, help='Path of samtools ', required=False, default='samtools')
parser.add_argument('-bl', '--blat_dir',type=str, help='Path of blat ', required=False, default='blat')
parser.add_argument('-ei', '--editing_index',type=str, help='Path of rna editing index tool ', required=False, default='RNAEditingIndex1.1')

EOF


SCRIPT_DIR="$(dirname "$(readlink -f "$0")")"

# run res scanner
${SCRIPT_DIR}/OneClick-RES-Scanner/easy_run_RES_scanner.sh -m $MAX_WORKERS -d $PROJECT_DIR -f $FASTQ_DIR -su $FASTQ_SUFFIX -dn $DNA_FILE -ln $DNA_LENGTH -lr $RNA_LENGTH -r $REF_FILE -a $ADDITIONAL_ARGUMENTS -b $BWA_DIR -s $SAMTOOLS_DIR -bl $BLAT_DIR

