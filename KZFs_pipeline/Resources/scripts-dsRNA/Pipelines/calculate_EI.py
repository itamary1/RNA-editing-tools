import pandas as pd
import sys
__author__ = "Itamar Twersky"

# input argumnts:
# sys.argv[1] = StrandDerivingCountsPerRegion.csv
# sys.argv[2] = True/False if to create total EI for all the regions tougether - if you have overlapping regions its not right to do that
# sys.argv[3] = deserved output file path for total  EI for all the regions tougether
'''
caclulate editing index given mismatch 
for StrandDerivingCountsPerRegion.csv table(should be an subtable of one sample)
MM format is X2X - for example "A2G"
'''
def calculateMM(MM, table):
    edited_count_col="NumOf"+MM+"Mismatches"
    match_count_col="NumOf"+MM.split('2')[0]
    edited_sum=table[edited_count_col].sum()
    return (edited_sum/(edited_sum+table[match_count_col].sum()))*100

    

# input argumnts:
# sys.argv[1] = input StrandDerivingCountsPerRegion.csv
# sys.argv[2] = True/False if to create total EI for all the regions tougether - if you have overlapping regions its not right to do that
# sys.argv[3] = deserved output file path for total  EI for all the regions tougether (can be empty if sys.argv[2] was set to false)
if __name__ == '__main__':
    input_table=pd.read_csv(sys.argv[1])
    mm_list=['A2C','A2G','A2T','C2A','C2G','C2T']
    # create editing index for every region in StrandDerivingCountsPerRegion
    for mm in mm_list:
        edited_count_col="NumOf"+mm+"Mismatches"
        match_count_col="NumOf"+mm.split('2')[0]
        input_table[mm+"EditingIndex"]= (input_table[edited_count_col]/(input_table[edited_count_col]+input_table[match_count_col]))*100
    input_table.to_csv(sys.argv[1],index=False)
    if(sys.argv[2]=='True'):
        assert len(sys.argv)==4, 'you sould give 3 argumnts if you want total EI'
        samples_list=list(input_table.Sample.unique())
        result_l=[]
        for sample in samples_list:
            sample_table=input_table[input_table['Sample'] == sample]
            first_line=sample_table.iloc[0]
            StrandDecidingMethod = first_line.StrandDecidingMethod
            Group = first_line.Group
            SamplePath = first_line.SamplePath
            Group=('Editing_Index_Unknown')
            index_dict={'StrandDecidingMethod' : StrandDecidingMethod, 'Group':Group, 'Sample':sample, 'SamplePath':SamplePath}
            for mm in mm_list:
                indexRes=calculateMM(mm,sample_table)
                index_dict[mm+"EditingIndex"] = "{:.14f}".format(indexRes)
            result_l.append(index_dict)
        final_table=pd.DataFrame(result_l)
        output_file=sys.argv[3]
        final_table.to_csv(output_file,index=False)
