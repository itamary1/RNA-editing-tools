import sys
import pandas as pd
import os
from distutils.dir_util import copy_tree
import subprocess
import shutil

nextflowSTE_config = "/home/alu/twerski/Scripts/Nextflow/Special_pipelines/Configs/Dockers/SubP_configs/salmonTE.nf.docker.config"
nextflowSTE_nf = "/home/alu/twerski/Scripts/Nextflow/Special_pipelines/Subpipelines/salmonTE.nf"

def move_out_files(condition_table,STE_quant_dir,out_dir):
    if not os.path.exists(out_dir):
        os.mkdir(out_dir)
    with open(condition_table, 'r') as fd:
       samples = [l.split(',')[0] for i,l in enumerate(fd.read().splitlines()) if i > 0]
    # sample_dirs = [os.path(STE_quant_dir,d) for d in os.listdir(STE_quant_dir) if os.path.isdir(d) and d in samples]
    [copy_tree(os.path.join(STE_quant_dir,sample), os.path.join(out_dir,sample))  for sample in samples]
    # mapp info
    mappInf=pd.read_csv(os.path.join(STE_quant_dir,'MAPPING_INFO.csv'))
    mappInf=mappInf[mappInf['SampleID'].isin(samples)]
    mappInf.to_csv(os.path.join(out_dir,'MAPPING_INFO.csv'),index=False)
    # expr
    expr=pd.read_csv(os.path.join(STE_quant_dir,'EXPR.csv'))
    samples.insert(0,'TE')
    expr=expr[samples]
    expr.to_csv(os.path.join(out_dir,'EXPR.csv'),index=False)
    # conditions
    cond=pd.read_csv(os.path.join(STE_quant_dir,'condition.csv'))
    cond=cond[cond['SampleID'].isin(samples)]
    cond.to_csv(os.path.join(out_dir,'condition.csv'),index=False)
    # clades
    shutil.copyfile(os.path.join(STE_quant_dir,'clades.csv'), os.path.join(out_dir,'clades.csv'))

def split_condition(condition_table,out_dir):
    cond_t_df = pd.read_csv(condition_table)
    conditions=cond_t_df.Condition.unique()
    couples = [(conditions[i],conditions[i+1]) for i in range(len(conditions)) if i+1 < len(conditions)]
    if not os.path.exists(out_dir):
        os.mkdir(out_dir)
    conditions_splited=[]
    for i,coup in enumerate(couples):
        coup_table = cond_t_df[cond_t_df.Condition.isin(coup)]
        outP=os.path.join(out_dir,"conditions{}.csv".format(i))
        coup_table.to_csv(outP,index=False)
        conditions_splited.append((outP,coup))
    return conditions_splited

def splitAndRun(condition_table,STE_quant_dir,out_dir):
    condition_subTabels_dir = os.path.join(out_dir,"condition_subTabels")
    if not os.path.exists(out_dir):
        os.mkdir(out_dir)
    os.mkdir(condition_subTabels_dir)
    condition_subTabels_list = split_condition(condition_table,condition_subTabels_dir)
    # condition_subTabels_list is list og tuple(two_conditions_list, condition_file)
    for cTable in condition_subTabels_list:
        # cTable is tuple(two_conditions_list, condition_file)
        subDir = os.path.join(out_dir,'_'.join(cTable[1]))
        os.mkdir(subDir)
        subQ_dir=os.path.join(subDir,'quant_files')
        move_out_files(cTable[0],STE_quant_dir,subQ_dir)
        os.chdir(subDir)
        p = subprocess.run(f"nextflow -bg -c {nextflowSTE_config} run {nextflowSTE_nf} --salmonTE_outdir {subDir} --test_step_only --control_condition {cTable[1][0]} --quant_dir {subQ_dir} --conditions_file {cTable[0]} &> {subDir}/nfRun.log.txt",shell=True)


if __name__ == "__main__":
    if(sys.argv[1] == "--help"):
        print("use with: subset condition_table STE_quant_dir out_dir")
        print("or use with: splitCond condition_table out_dir")
    if(sys.argv[1] == "splitCond"):
        condition_table = sys.argv[2]
        out_dir = sys.argv[3]
        split_condition(condition_table,out_dir)
    if(sys.argv[1] == "subset"):
        condition_table = sys.argv[2]
        STE_quant_dir = sys.argv[3]
        out_dir = sys.argv[4]
        move_out_files(condition_table,STE_quant_dir,out_dir)
    if(sys.argv[1] == "splitAndRun"):
        condition_table = sys.argv[2]
        STE_quant_dir = sys.argv[3]
        out_dir = sys.argv[4]
        splitAndRun(condition_table,STE_quant_dir,out_dir)



   
    