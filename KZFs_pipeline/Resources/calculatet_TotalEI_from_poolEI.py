
import sys
import pandas as pd
import os

'''
arguments:
    pool_final_tabel = sys.argv[1]
    num_gene_to_take = int(sys.argv[2])
    out_file = sys.argv[3]

'''


mm_list=['A2C','A2G','A2T','C2A','C2G','C2T']

def calculateMM(MM, table):
    ref,mm=MM.split('2')
    edited_count_col=mm+'_ReferenceBase'+ref
    match_count_col=ref+'_ReferenceBase'+ref
    edited_sum=table[edited_count_col].sum()
    return (edited_sum/(edited_sum+table[match_count_col].sum()))*100

def calculateEI_gene_list(pool_final_tabel,gene_list,filter_cov,debugC=False):
    pool_final_tabel = pool_final_tabel[pool_final_tabel.AGCoveragePerSite > filter_cov]
    genes_tabel = pool_final_tabel[pool_final_tabel.Region.isin(gene_list)]
    if debugC:
        genes_tabel.to_csv(f"/private10/Projects/KZNF_itamar/KZFs_pipeline/Runs/GSE190548_6embryoStages/Results/Editing_index/KZF_CDS_pooledEI/top{len(gene_list)}_genes_summaryTable.csv",index=False)
    samples_list=list(genes_tabel.Sample.unique())
    result_l=[]
    for sample in samples_list:
        sample_table=genes_tabel[genes_tabel['Sample'] == sample]
        StrandDecidingMethod = "bed_file"
        Group = "noGroup"
        SamplePath = "unknown"
        Group=('Editing_Index_Unknown')
        index_dict={'StrandDecidingMethod' : StrandDecidingMethod, 'Group':Group, 'Sample':sample, 'SamplePath':SamplePath}
        for mm in mm_list:
            indexRes=calculateMM(mm,sample_table)
            index_dict[mm+"EditingIndex"] = "{:.14f}".format(indexRes)
        result_l.append(index_dict)
    final_table=pd.DataFrame(result_l)
    return final_table

def get_top_genes(top_genes_table,num_gene_to_take):
    top_genes =(top_genes_table['Region'].unique())[:num_gene_to_take]
    print(len(top_genes),"got", num_gene_to_take)
    return top_genes


if __name__ == "__main__":
    pool_final_tabel = sys.argv[1]
    num_gene_to_take = int(sys.argv[2])
    top_genes_table = sys.argv[3]
    out_file = sys.argv[4]
    pool_final_tabelDf=pd.read_csv(pool_final_tabel)
    top_genes_tableDf = pd.read_csv(top_genes_table)
    top_genes=get_top_genes(top_genes_tableDf,num_gene_to_take)
    final_table = calculateEI_gene_list(pool_final_tabelDf,top_genes,500,False)
    final_table.to_csv(out_file,index=False)

