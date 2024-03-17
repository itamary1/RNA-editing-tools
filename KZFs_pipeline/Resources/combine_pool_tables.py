import sys
import csv
import os
import pandas as pd

# this script will create final table of per gene pooling scripts
# it will select its strand by the bed6 file
# pooling_result_dir = sys.argv[1]
# bed6_file = sys.argv[2]
# out_table_path = sys.argv[3]

mm_list=['A2C','A2G','A2T','C2A','C2G','C2T']

def combine_strand_tables(pooling_result_dir,type):
    all_tables = []
    for table in [t for t in os.listdir(pooling_result_dir) if t.startswith('EditedRegionsSites.pooled') and t.endswith(type+'.csv')]:
        sample_name = str(table).removeprefix('EditedRegionsSites.pooled.').removesuffix(f".{type}.csv")
        tableP=os.path.join(pooling_result_dir,table)
        sample_df = pd.read_csv(tableP)
        sample_df['Sample'] = sample_name
        EI_list =  [mm+'EditingIndexNoSNP' for mm in mm_list]
        count_match_l = ['C_ReferenceBaseC','A_ReferenceBaseA']
        mm_count_l = [mm.split('2')[1] +'_ReferenceBase' +mm.split('2')[0] for mm in mm_list]
        selectedColumns = ['Sample','Region','AGCoveragePerSite']
        selectedColumns.extend(EI_list)
        selectedColumns.extend(count_match_l)
        selectedColumns.extend(mm_count_l)
        sample_df=sample_df[selectedColumns]
        all_tables.append(sample_df)
    return pd.concat(all_tables)



    

if __name__ == "__main__":
    pooling_result_dir = sys.argv[1]
    bed6_file = sys.argv[2]
    out_table_path = sys.argv[3]
    AG_tables = combine_strand_tables(pooling_result_dir,'A2G')
    TC_tables = combine_strand_tables(pooling_result_dir,'T2C')
    all_genes =[]
    with open(bed6_file, 'r') as fd:
        for gene in csv.DictReader(fd, fieldnames=['Chr','Start','End','GeneName','bla','Strand'],delimiter="\t"):
            assert gene['Strand'] in ['-','+'], 'strand nor -/+'
            strandTables = AG_tables if gene['Strand'] == '+' else TC_tables
            all_genes.append(strandTables[strandTables.Region == gene['GeneName']])
    final_table =pd.concat(all_genes)
    namesDict={name:str(name).removesuffix('NoSNP') for name in final_table.columns}
    final_table.rename(columns=namesDict,inplace=True)
    final_table.to_csv(out_table_path,index=False)
