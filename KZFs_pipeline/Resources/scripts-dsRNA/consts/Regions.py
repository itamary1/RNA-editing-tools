import os
if __name__ == '__main__':
    sys.path.append(os.path.dirname(os.path.dirname(__file__)))
from common.files_general import read_json

# region Consts
consts_json = read_json(os.path.join(os.path.dirname(__file__), "Regions2.json"))
# get
SPLIT_BY_VAR=consts_json['SPLIT_BY_VAR']
COLNAMES = consts_json['COLNAMES']
GROUPING_VARS = consts_json['GROUPING_VARS']
SPLIT_OUT_PREFIX= consts_json['SPLIT_OUT_PREFIX']
SPLIT_OUT_SUFFIX= consts_json['SPLIT_OUT_SUFFIX']
POOL_OUT_PREFIX= consts_json['POOL_OUT_PREFIX']
POOL_OUT_SUFFIX= consts_json['POOL_OUT_SUFFIX']
# endregion

# region Scripts and commands
# paths relative to current file
CHANGE_FILE_NAMES_CMD="find %(input_dir)s -maxdepth 2 -type f -name StrandDerivingCountsPerRegion.csv -exec bash -c \'mv $0 ${0/StrandDerivingCountsPerRegion.csv/StrandDerivingCountsPerRegion_notPerSample.csv}\' {} \;"
# split file by chromosomes
SPLIT_BY_CHR_SCRIPT=os.path.join(os.path.dirname(os.path.dirname(__file__)), "Processing", "Editing", "split_StrandDerivingCountsPerRegion_file_by_chromosomes.R")
SPLIT_BY_CHR_CMD="%(r)s %(split_by_chr)s " + "-i %(input_dir)s -o %(input_dir)s" + (" -is StrandDerivingCountsPerRegion.csv -p %(prefix)s -os %(suffix)s -f %(split_var)s --wanted_columns %(wanted_cols)s" % {'prefix' : SPLIT_OUT_PREFIX, 'suffix' : SPLIT_OUT_SUFFIX, 'split_var' : SPLIT_BY_VAR, 'wanted_cols' : COLNAMES})
# pool splitted files by chromosomes
POOL_BY_CHR_SCRIPT=os.path.join(os.path.dirname(os.path.dirname(__file__)), "Processing", "Editing", "Region_perRegionPerSample_GTExSubset_preprocess.withSignal.v3.R")
POOL_BY_CHR_CMD="%(r)s %(pool_by_chr)s -i %(input_dir)s -o %(output_dir)s -g %(group_file)s" + " -c %(wanted_cols)s -v %(group_vars)s -p %(input_prefix)s -s %(input_suffix)s -op %(output_prefix)s -os %(output_suffix)s" %{'wanted_cols' : SPLIT_BY_VAR + " " + COLNAMES, 'group_vars' : GROUPING_VARS, 'input_prefix' : SPLIT_OUT_PREFIX, 'input_suffix' : SPLIT_OUT_SUFFIX, 'output_prefix' : POOL_OUT_PREFIX, 'output_suffix' : POOL_OUT_SUFFIX}
# awk commands
COMBINE_FILES_CMD="awk '{if(NR==1 || (FNR > 1)) print $0}' %(output_dir)s/stats/" + "%(output_prefix)s.chr*.%(output_suffix)s" % {'output_prefix' : POOL_OUT_PREFIX, 'output_suffix' : POOL_OUT_SUFFIX} + "> %(output_dir)s/stats/" + "%(output_prefix)s.%(output_suffix)s" % {'output_prefix' : POOL_OUT_PREFIX, 'output_suffix' : POOL_OUT_SUFFIX}
COMBINE_COUNTS_CMD='awk \'BEGIN{FS=OFS=","; print "AluElement,GroupCount_Edited,GroupCount_NotEdited,SampleCount_Edited,SampleCount_NotEdited,Groups_Edited,Groups_NotEdited"} {if(FNR > 1 && NF == 7) {print $0} else if (FNR > 1 && NF < 7) {print $1, $2, "NA", $3, "NA", $4, "NA"}}\' %(output_dir)s/stats/CountDetails' + "%(output_prefix)s.chr*.%(output_suffix)s" % {'output_prefix' : POOL_OUT_PREFIX, 'output_suffix' : POOL_OUT_SUFFIX} + ' >  %(output_dir)s/stats/CountDetails' + "%(output_prefix)s.%(output_suffix)s" % {'output_prefix' : POOL_OUT_PREFIX, 'output_suffix' : POOL_OUT_SUFFIX}
# endregion