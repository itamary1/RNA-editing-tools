#!/usr/bin/env bash

#for ruunnig All RES parts dou to configs files in dir

ARGPARSE_DESCRIPTION="ruunnig All RES parts dou to configs files in dir - this script assunmmig you have 1 dna fastq file and many rna fastq files"      
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
EOF

SCRIPT_DIR="$(dirname "$(readlink -f "$0")")"

# create config file
mkdir -p $PROJECT_DIR/Configs
${SCRIPT_DIR}/Resources/create_res_configs.sh -o $PROJECT_DIR/Configs -f ${FASTQ_DIR} -su $FASTQ_SUFFIX  -dn $DNA_FILE -ln $DNA_LENGTH -lr $RNA_LENGTH

counter=1;

ls $CONFIGS_DIR | sort -n |while read -r REL_CONF_FILE
do
	CONF_FILE=${CONFIGS_DIR}${REL_CONF_FILE}
	SAMPLE_NAME="${REL_CONF_FILE%\_*}"
	OUT_DIR=${PROJECT_DIR}"/"${SAMPLE_NAME}"/"
	mkdir -p $OUT_DIR
	cd $OUT_DIR
	nohup "${SCRIPT_DIR}/Resources/run_RES_oneConf_allParts.sh" -o $OUT_DIR -c $CONF_FILE -r $REF_FILE -a "$ADDITIONAL_ARGUMENTS" -b $BWA_DIR -s $SAMTOOLS_DIR -l $BLAT_DIR &
	pids[${counter}]=$!
	if (($counter%MAX_WORKERS == 0))
	then
		echo "waiting"
		# wait for all pids
		for pid in ${pids[*]}; do
			wait $pid
			echo "sample proccess pid: $pid Exit status: $?"
		done
		pids=()
	fi
	counter=$(($counter+1))

done

# wait for all 
wait

echo "run reas scanner on all samples done"
