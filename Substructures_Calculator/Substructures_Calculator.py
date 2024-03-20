import pandas as pd
import sys
import RNAstracture
import bedtoolsIT
import traceback

def print_help_mesage():
    massage='''

    ----------------------------------------------------------------------------------------------------------

                                                    Substructures Calculator

    -----------------------------------------------------------------------------------------------------------
    This tool are for getting Substructures information of editing sites

    Its accept csv table of editing sites with columns of 'Coord' and 'Chr'
    
    the input arguments are positional arguments - make sure its the correct order

    inputs argumnts:
    genome - genome in fasta format
    editing_sites_table  - csv table of editing sites with columns of 'Coord' and 'Chr'
    data_tabels - FOLD datatables
    FOLD - path to FOLD program
    bpRNA - path to bprna program
    out_table_path - path to write output table
    window (optional)- window bp in each side of every editing site to take for calculating all the substructures around the editing site - defualt 400

    an example runing command:
    
    python Substructures_Calculator.py $genome $editing_sites_table $data_tabels $FOLD $bpRNA $out_table_path $window

    See full runing example at usage_example.sh



    previous dependecy - you need to have:
        -python 3.9 with pandas and SeqIO
        -Fold from RNAstructure package and its datatables https://rna.urmc.rochester.edu/RNAstructureDownload.html
        -perl with Graph.pm
        -bpRNA https://github.com/hendrixlab/bpRNA


    '''
    print(massage)



def get_seq(site, window,fasta_getter):
    """will return sequnce around site at site window in every side

    Args:
        site (a one site from df): site from df.iterrows()
        window (int): size of window

    Returns:
        str: sequnce around site
    """
    start = int(site["Coord"]) - window - 1
    end = int(site["Coord"]) + window - 1
    seq = fasta_getter.get_fasta(site["Chr"], start, end)
    return seq

if __name__ =="__main__":
    if len(sys.argv) < 6 or sys.argv[1]=='-h' or sys.argv[1]=='--help':
        print_help_mesage()
        exit(0)
    try:
        genome_file=sys.argv[1]
        in_path=sys.argv[2]
        data_tables=sys.argv[3]
        Fold=sys.argv[4]
        bpRNA=sys.argv[5]
        out_path=sys.argv[6]
        window= 400 if len(sys.argv) < 7 else sys.argv[7]
        # create an object for getting fasta sequnces from a genome base on coordinates
        # the objext will read the whole genome to the memory for fast respone to requests
        fasta_getter = bedtoolsIT.genome_reader(genome_file)
        # read editing sites table
        editing_sites_table=pd.read_csv(in_path)
        # get sequece of 400bp in each side of every editing site
        seq_list = [get_seq(site,window,fasta_getter) for i,site in editing_sites_table.iterrows()]
        # create Near_site_folding_calculator
        nc = RNAstracture.Near_site_folding_calculator(data_tables,Fold,bpRNA)
        # calculate substructurs for every editing site
        fold_results = nc.calculate_substrcture_seq_list(seq_list)
        # add the folding result to the editing sites table
        editing_sites_table["big_wind_seq"] = seq_list
        editing_sites_table["big_wind_f_energy"] = [res["big_f_e"] if res["big_f_e"] else np.nan for res in fold_results]
        editing_sites_table["big_wind_struct"] = [res["big_s"] if res["big_s"] else np.nan for res in fold_results]
        editing_sites_table["small_wind_seq"] = [res["small_seq"] if res["small_seq"] else np.nan for res in fold_results]
        editing_sites_table["small_seq_coords"] = [
            res["small_seq_coords"] if res["small_seq_coords"] else np.nan for res in fold_results
        ]
        editing_sites_table["small_wind_f_energy"] = [res["small_f_e"] if res["small_f_e"] else np.nan for res in fold_results]
        editing_sites_table["small_wind_struct"] = [res["small_s"] if res["small_s"] else np.nan for res in fold_results]
        editing_sites_table.to_csv(out_path,index=False)
    except Exception as e:
        print("error occur during runing")
        print("make sure you download Fold and bpRNA and set pathes correctlys")
        print("your error:")
        print(traceback.format_exc())
        print_help_mesage()

