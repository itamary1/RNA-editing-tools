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
parser.add_argument('-b', '--bwa_path',type=str, help='Path of BWA ', required=False, default='bwa')
parser.add_argument('-s', '--samtools_path',type=str, help='Path of samtools ', required=False, default='samtools')
parser.add_argument('-bl', '--blat_path',type=str, help='Path of blat ', required=False, default='blat')
parser.add_argument('-bt', '--bedtools_path',type=str, help='Path of bedtools ', required=False, default='bedtools')
parser.add_argument('-py', '--python_path',type=str, help='Path of python with pandas and numpy ', required=False, default='python')
parser.add_argument('-rq', '--refseq_file', type=str, help='Path of refseq gene expression file', required=True)
parser.add_argument('-ei', '--editing_index',type=str, help='Path of rna editing index tool ', required=False, default='RNAEditingIndex')

EOF


SCRIPT_DIR="$(dirname "$(readlink -f "$0")")"

# run res scanner
${SCRIPT_DIR}/OneClick-RES-Scanner/easy_run_RES_scanner.sh -m $MAX_WORKERS -d $PROJECT_DIR -f $FASTQ_DIR -su $FASTQ_SUFFIX -dn $DNA_FILE -ln $DNA_LENGTH -lr $RNA_LENGTH -r $REF_FILE -a $ADDITIONAL_ARGUMENTS -b $BWA_PATH -s $SAMTOOLS_PATH -bl $BLAT_PATH

reagions_dir=$PROJECT_DIR"/All_sites"
tables_dir=${reagions_dir}"/final_tabels"

# a file for all site positions - clear it if exist
regions=${reagions_dir}"/All_sites.bed"
mkdir -p $tables_dir
truncate -s 0 $regions
# find all res results files, parse them, write final tabels to a dir and create bed file from all the genomics positions
for res_resu in $(find $PROJECT_DIR -name "RES_final_result.txt" -type f); do
    file_name=${res_resu%/identification/*}
    file_name=${file_name##*/}
    out_file_path=$tables_dir"/"${file_name}"_final_result.tsv"
    $PYTHON_PATH $SCRIPT_DIR"/"parse_res_table.py $res_resu $(wc -m < $REF_FILE) 0.05 $out_file_path
    cat $out_file_path | awk 'BEGIN{FS="\t";OFS="\t"};NR>1{print $1,$2-1,$2}' >> $regions
done
regions_merged_sorted=${PROJECT_DIR}"/All_sites/All_sites.sorted.merged.bed"
$BEDTOOLS_PATH sort -i $regions | bedtools merge > $regions_merged_sorted
# get all bams
all_RNA_bams_dir=$PROJECT_DIR"/RNA_bams_links"
mkdir -p $all_RNA_bams_dir
find $PROJECT_DIR -name "RNA" -type d | xargs -I {} find {} -name "*bwa.bam" -type f | xargs -I BF ln -s BF ${all_RNA_bams_dir}/


editing_index_dir=$PROJECT_DIR"/RNA_editing_index"
mkdir -p $editing_index_dir
$EDITING_INDEX -d ${all_RNA_bams_dir} -l ${editing_index_dir}/LOGS/ -o ${editing_index_dir}/Output_pileup/ -os "${editing_index_dir}/Index_output/" --genome UserProvided -gf $REF_FILE --refseq $REFSEQ_FILE -rb $regions_merged_sorted --bam_files_suffix ".bwa.bam" --follow_links