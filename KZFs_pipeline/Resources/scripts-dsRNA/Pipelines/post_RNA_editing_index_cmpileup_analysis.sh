# create argsparse
ARGPARSE_DESCRIPTION="Pool post-processing for RNA editing index"      # this is optional
source /private/common/Software/BashLibs/argparse-bash/argparse.bash || exit 1
argparse "$@" <<EOF || exit 1
parser.add_argument('-i', '--input_dir', type=str, help='Path of RNA editing index directory (where cmpileup directory resides)', required=True)
parser.add_argument('-o', '--output_dir',type=str, help='Path of pooled results directory', required=True)
parser.add_argument('--group_file',nargs='+', type=str, help='CSV with sample to group. Will pool samples according to given groups, sample count per group will be calculated according to this file. Group must include fragment of file path (file name or directories). Format: Sample,Group',required=True)
parser.add_argument('-s', '--suff_cmpileup',type=str, help='cmpileup suffix', default="_ucscHg38Alu.bed.gz_mpileup.cmpileup")
parser.add_argument('-e', '--email',type=str, help='Email of user', default="ronif10+levanonlab@gmail")
EOF

# stop if failed
set -e
trap 'last_command=$current_command; current_command=$BASH_COMMAND' DEBUG
trap 'echo ""${last_command}" command filed with exit code $?."' EXIT

###############################
####    CONSTANT SUFFIXES   ###
###############################
python_cmd=python3.6
r_cmd=Rscript
sum_cmpileup_allMM_script_cmd=/home/alu/twerski/Scripts/Nextflow/Special_pipelines/Resources/scripts-dsRNAProcessing/Cmpileup/sum_cmpileup_allMM.py
split_summed_cmpileup_by_SNP=/home/alu/twerski/Scripts/Nextflow/Special_pipelines/Resources/scripts-dsRNAProcessing/Cmpileup/summed_cmpileup_filterSNP.sh
pool_info_script=/home/alu/twerski/Scripts/Nextflow/Special_pipelines/Resources/scripts-dsRNAProcessing/Cmpileup/all_ES_in_regions_analysis.R

###############################
####    FLATTEN CMPILEUPS   ###
###############################
echo "DEBUG: Pooling cmpileup by group"
# "flatten" cmpileup across groups by summing each group
for g in $(awk 'NR>1' $GROUP_FILE | cut -f 2 -d , | sort -u | tr "\n" " "); do
	echo "DEBUG: Pooling cmpileup by group - ${g}"
	echo "DEBUG: Pooling cmpileup by group - ${g}: $python_cmd $sum_cmpileup_allMM_script_cmd -d $INPUT_DIR -s $SUFF_CMPILEUP -o $INPUT_DIR/Combined.$g.$SUFF_CMPILEUP.pooled.csv --sample_list $(grep $g $GROUP_FILE | cut -f 1 -d , | sort -u | tr "\n" " ")"
	$python_cmd $sum_cmpileup_allMM_script_cmd -d $INPUT_DIR -s $SUFF_CMPILEUP -o $INPUT_DIR/Combined.$g.$SUFF_CMPILEUP.pooled.csv --sample_list $(grep $g $GROUP_FILE | cut -f 1 -d , | sort -u | tr "\n" " ") > post_RNA_editing_index_cmpileup_analysis.pool_cmpileup.$g.log
done

#################################
####    SPLIT SNP POSITIONS   ###
#################################
echo "DEBUG: Splitting into SNP and non-SNP sites"
### split into SNP and non-SNP points for each file, strand-corrected #############################################
# per group
for g in $(awk 'NR>1' $GROUP_FILE | cut -f 2 -d , | sort -u | tr "\n" " "); do
	echo "DEBUG: Splitting into SNP and non-SNP sites - ${g}"
	echo "DEBUG: Splitting into SNP and non-SNP sites - ${g}: sh $split_summed_cmpileup_by_SNP -i $INPUT_DIR/Combined.$g.$SUFF_CMPILEUP.pooled.csv -o $INPUT_DIR/Combined.$g.$SUFF_CMPILEUP.pooled.noSNP.csv -os $INPUT_DIR/Combined.$g.$SUFF_CMPILEUP.pooled.onlySNP.csv"
	sh $split_summed_cmpileup_by_SNP -i $INPUT_DIR/Combined.$g.$SUFF_CMPILEUP.pooled.csv -o $INPUT_DIR/Combined.$g.$SUFF_CMPILEUP.pooled.noSNP.csv -os $INPUT_DIR/Combined.$g.$SUFF_CMPILEUP.pooled.onlySNP.csv > post_RNA_editing_index_cmpileup_analysis.split_summed_cmpileup_by_SNP.$g.log
done

##############################
####    POOL INFORMATION   ###
##############################
echo "DEBUG: Pooling editing site information"
echo "DEBUG: Pooling editing site information: $r_cmd $pool_info_script -i $INPUT_DIR -o $OUTPUT_DIR -g $GROUP_FILE"
$r_cmd $pool_info_script -i $INPUT_DIR -o $OUTPUT_DIR -g $GROUP_FILE > post_RNA_editing_index_cmpileup_analysis.pool_info.log

echo "INFO: Done"