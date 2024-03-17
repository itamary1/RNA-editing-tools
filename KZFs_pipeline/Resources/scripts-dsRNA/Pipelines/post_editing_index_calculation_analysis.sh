# create argsparse
ARGPARSE_DESCRIPTION="Pool post-processing for RNA editing index"      # this is optional
source /private/common/Software/BashLibs/argparse-bash/argparse.bash || exit 1
argparse "$@" <<EOF || exit 1
parser.add_argument('-i', '--input_dir', type=str, help='Path of RNA editing index directory (where cmpileup directory resides)', required=True)
parser.add_argument('-o', '--output_dir',type=str, help='Path of pooled results directory', required=True)
parser.add_argument('--group_file',nargs='+', type=str, help='CSV with sample to group. Will pool samples according to given groups, sample count per group will be calculated according to this file. Format: Sample,Group',required=True)
parser.add_argument('-e', '--email',type=str, help='Email of user', default="ronif10+levanonlab@gmail")
EOF

# stop if failed
set -e
trap 'last_command=$current_command; current_command=$BASH_COMMAND' DEBUG
trap 'echo ""${last_command}" command filed with exit code $?."' EXIT


###############################
####    CONSTANT SUFFIXES   ###
###############################
r_cmd=Rscript
split_by_chr_script=/home/alu/twerski/Scripts/Nextflow/Special_pipelines/Resources/scripts-dsRNAProcessing/Editing/split_StrandDerivingCountsPerRegion_file_by_chromosomes.R
pool_by_chr_script=/home/alu/twerski/Scripts/Nextflow/Special_pipelines/Resources/scripts-dsRNAProcessing/Editing/Region_perRegionPerSample_GTExSubset_preprocess.withSignal.R

#################################
####   SPLIT REGIONS BY CHR   ###
#################################
echo "DEBUG: changing top-level StrandDerivingCountsPerRegion.csv file names"
# First change names of files that are grouped per-tissue (as they will be caught by the suffix)
find $INPUT_DIR -maxdepth 2 -type f -name StrandDerivingCountsPerRegion.csv -exec bash -c 'mv $0 ${0/StrandDerivingCountsPerRegion.csv/StrandDerivingCountsPerRegion_notPerSample.csv}' {} \;
echo "Names of top-level StrandDerivingCountsPerRegion.csv changed"

# run uniting script
echo "DEBUG: splitting information by chromosomes"
$r_cmd $split_by_chr_script -i $INPUT_DIR -is StrandDerivingCountsPerRegion.csv -o $INPUT_DIR -p StrandDerivingCountsPerRegion -os preprocessed.csv -f GenomicRegion --wanted_columns Sample AluElement Length SenseGeneCommonName SenseGeneRefSeqID AntisenseGenomicPosition SenseGenomicPosition AntisenseGeneRefSeqID AntisenseGeneCommonName A2GEditingIndex C2TEditingIndex TotalCoverageAtAllPositions MeanCoveragePerRegion IndexedMismatchesOfA2G IndexedCanonicalOfA2G IndexedMismatchesOfC2T IndexedCanonicalOfC2T NumOfIndexedMismatchesSitesOfA2G NumOfIndexedOverallSitesOfA2G NumOfIndexedMismatchesSitesOfC2T NumOfIndexedOverallSitesOfC2T > post_editing_index_calculation_analysis.split_by_chr.log

echo "DEBUG: splitting information by chromosomes for other MM"
# run uniting script other MM - part 1
$r_cmd $split_by_chr_script -i $INPUT_DIR -is StrandDerivingCountsPerRegion.csv -o $INPUT_DIR -p StrandDerivingCountsPerRegion -os preprocessed.otherMM1.csv -f GenomicRegion --wanted_columns Sample AluElement Length TotalCoverageAtAllPositions IndexedCanonicalOfC2A IndexedCanonicalOfG2C IndexedCanonicalOfA2G IndexedMismatchesOfC2G IndexedMismatchesOfC2A IndexedCanonicalOfG2T IndexedMismatchesOfG2C IndexedMismatchesOfG2T NumOfIndexedOverallSitesOfC2G NumOfIndexedOverallSitesOfC2A IndexedCanonicalOfC2G NumOfIndexedOverallSitesOfG2C NumOfIndexedOverallSitesOfG2T NumOfIndexedOverallSitesOfA2G NumOfIndexedMismatchesSitesOfG2C NumOfIndexedMismatchesSitesOfG2T IndexedMismatchesOfA2G NumOfIndexedMismatchesSitesOfC2A NumOfIndexedMismatchesSitesOfC2G NumOfIndexedMismatchesSitesOfA2G IndexedCanonicalOfG2A NumOfIndexedOverallSitesOfG2A IndexedMismatchesOfG2A NumOfIndexedMismatchesSitesOfG2A > post_editing_index_calculation_analysis.split_by_chr.part1.log

# run uniting script other MM - part 2
$r_cmd $split_by_chr_script -i $INPUT_DIR -is StrandDerivingCountsPerRegion.csv -o $INPUT_DIR -p StrandDerivingCountsPerRegion -os preprocessed.otherMM2.csv -f GenomicRegion --wanted_columns Sample AluElement Length TotalCoverageAtAllPositions IndexedCanonicalOfA2G IndexedCanonicalOfA2C IndexedCanonicalOfA2T NumOfIndexedMismatchesSitesOfT2A NumOfIndexedMismatchesSitesOfT2G IndexedMismatchesOfA2C IndexedMismatchesOfA2T IndexedCanonicalOfT2C IndexedCanonicalOfT2A NumOfIndexedMismatchesSitesOfT2C NumOfIndexedOverallSitesOfA2G NumOfIndexedOverallSitesOfA2C NumOfIndexedOverallSitesOfA2T NumOfIndexedOverallSitesOfT2G NumOfIndexedOverallSitesOfT2C NumOfIndexedOverallSitesOfT2A IndexedMismatchesOfA2G IndexedMismatchesOfT2G IndexedMismatchesOfT2A IndexedMismatchesOfT2C NumOfIndexedMismatchesSitesOfA2T NumOfIndexedMismatchesSitesOfA2C NumOfIndexedMismatchesSitesOfA2G IndexedCanonicalOfT2G > post_editing_index_calculation_analysis.split_by_chr.part2.log

#################################
####     POOL INFORMATION    ###
#################################
echo "DEBUG: pooling editing information by chromosomes"
echo "DEBUG: pooling editing information by chromosomes: $r_cmd $pool_by_chr_script -i $INPUT_DIR -o $OUTPUT_DIR -g $GROUP_FILE"
# pool each chromosome and index
$r_cmd $pool_by_chr_script -i $INPUT_DIR -o $OUTPUT_DIR -g $GROUP_FILE > post_editing_index_calculation_analysis.pool_by_chr.log

#################################
####     COMBINE ALL FILES    ###
#################################
echo "DEBUG: combining editing info files to dataset-level"
# combine processed files:
# pooling information by region
awk '{if(NR==1 || (FNR > 1)) print $0}' $OUTPUT_DIR/stats/EditingPerRegionPerSample.chr*.pooled.withSignal.csv > $OUTPUT_DIR/stats/EditingPerRegionPerSample.pooled.withSignal.csv
for part in otherMM1 otherMM2; do
	# by region
	awk '{if(NR==1 || (FNR > 1)) print $0}' $OUTPUT_DIR/stats/EditingPerRegionPerSample.chr*.pooled.$part.csv > $OUTPUT_DIR/stats/EditingPerRegionPerSample.pooled.$part.csv
done

# consider differently, as some alt chromosomes only have 4 fields due to regions that were only expressed and edited
# (no "_NotEdited" columns)
# [result will be NA for no expression information and empty/0 for no editing information]
# by region
echo "DEBUG: Combining count detail information"
awk 'BEGIN{FS=OFS=","; print "AluElement,GroupCount_Edited,GroupCount_NotEdited,SampleCount_Edited,SampleCount_NotEdited,Groups_Edited,Groups_NotEdited"} {if(FNR > 1 && NF == 7) {print $0} else if (FNR > 1 && NF < 7) {print $1, $2, "NA", $3, "NA", $4, "NA"}}' $OUTPUT_DIR/stats/CountDetailsEditingPerRegionPerSample.chr*.pooled.withSignal.csv > $OUTPUT_DIR/stats/CountDetailsEditingPerRegionPerSample.pooled.withSignal.csv
echo "Combined processed count information"

echo "DEBUG: DONE"
# send final email
mail -s 'Pooling of ${INPUT_DIR} Has Finished' $email <<< 'Analysis for ${INPUT_DIR} into ${OUTPUT_DIR}/stats is done. Enjoy :)'
