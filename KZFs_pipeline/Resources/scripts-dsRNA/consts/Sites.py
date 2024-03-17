import os

# region Scripts and commands
POOLED_OUT_FORMAT = "%(output_dir)s/Combined.%(group)s.%(cmpileup)s.pooled.csv"
POOLED_OUT_NO_SNP_FORMAT = "%(output_dir)s/Combined.%(group)s.%(cmpileup)s.pooled.noSNP.csv"
POOLED_OUT_W_SNP_FORMAT = "%(output_dir)s/Combined.%(group)s.%(cmpileup)s.pooled.onlySNP.csv"
SUM_CMPILEUP_SCRIPT = os.path.join(os.path.dirname(os.path.dirname(__file__)), "Processing", "Cmpileup", "sum_cmpileup_allMM.py")
# SUM_CMPILEUP_CMD = '%(python)s %(sum_cmpileup_script)s -d %(input_dir)s -s %(cmpileup)s --sample_list %(samples)s -o ' + POOLED_OUT_FORMAT
# FORCE_SUM_CMPILEUP_CMD = SUM_CMPILEUP_CMD + " --force"
SPLIT_CMPILEUP_BY_SNP_SCRIPT=os.path.join(os.path.dirname(os.path.dirname(__file__)), "Processing", "Cmpileup", "summed_cmpileup_filterSNP.sh")
SPLIT_CMPILEUP_BY_SNP_CMD='%(shell)s %(split_summed_cmpileup_by_SNP)s ' + '-i %(input)s -o %(out_no_snp)s -os %(out_w_snp)s' % {'input' : POOLED_OUT_FORMAT, 
                                                                                                                                'out_no_snp' : POOLED_OUT_NO_SNP_FORMAT, 
                                                                                                                                'out_w_snp' : POOLED_OUT_W_SNP_FORMAT} + " --snps %(snps)s"
POOL_ES_SCRIPT=os.path.join(os.path.dirname(os.path.dirname(__file__)), "Processing", "Cmpileup", "all_ES_in_regions_analysis.R")
# POOL_ES_CMD="%(r)s %(pool_info_script)s -i %(input_dir)s -o %(output_dir)s -g %(group_file)s"
POOL_ES_CMD="%(r)s %(pool_info_script)s -ins %(input_no_snps)s -ios %(input_only_snps)s -o %(output_dir)s -g %(group_name)s -n %(group_size)s"
POOL_ES_FORCE_PARAM = " --force"
POOL_ES_ALTERNATIVE_REGION_ID = " --region_id_file %s"
# endregion