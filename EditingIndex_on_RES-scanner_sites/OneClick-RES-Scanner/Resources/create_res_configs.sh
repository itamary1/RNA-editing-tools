#!/usr/bin/env bash

#for creating RES's configs. one config for every fastq

#for Ffile in $(ls *fastq); do mv /private6/Projects/Yeast_Itamar_03_2022/FASTQ_Generation_2022-03-14_22_58_33Z-540314775/Unzipped_merged_fastq/${Ffile} $(echo ${Ffile#*S} | sed 's/_/S_/'); done

ARGPARSE_DESCRIPTION="ruunnig RES's all parts for one config"      # this is optional
source ~/argparse-bash/argparse.bash || exit 1
argparse "$@" <<EOF || exit 1
parser.add_argument('-o', '--out_dir', type=str, help='Full path of existing output directory', required=True)
parser.add_argument('-f', '--fastq_dir', type=str, help='Path of fastq dir', required=True)
parser.add_argument('-su', '--fastq_suffix',type=str, help='fastq files suffix', required=False, default='fastq')
parser.add_argument('-dn', '--dna_file', type=str, help='full Path of dna file', required=True)
parser.add_argument('-ln', '--dna_length',type=str, help='length of reads of DNA', required=False, default='82')
parser.add_argument('-lr', '--rna_length',type=str, help='length of reads of RNA', required=False, default='81')

EOF

cd $FASTQ_DIR

for Ffile in $(ls *${fastq_suffix} | sort -n)
do
	FULL_PATH_FASTQ=${FASTQ_DIR}/${Ffile}
	SAMPLE_NAME="${Ffile%%\.*}"
	NEW_FILE=${OUT_DIR}/${SAMPLE_NAME}_config.txt
	
	echo -e "DNA""\t""$SAMPLE_NAME""\t""s0""\t""$DNA_LENGTH""\t""$DNA_FILE" > $NEW_FILE
	echo -e "RNA""\t""$SAMPLE_NAME""\t""$SAMPLE_NAME""\t""$RNA_LENGTH""\t""$FULL_PATH_FASTQ" >> $NEW_FILE
	
done
