import pandas as pd
import sys
import RNAstracture
import bedtoolsIT

# directory of our current script
script_dir=sys.argv[0].removesuffix("usage_example.py")

# create an object for getting fasta sequnces from a genome base on coordinates
# the objext will read the whole genome to the memory for fast respone to requests
fasta=script_dir+"first10000_sacCer3.fa"
fasta_getter = bedtoolsIT.genome_reader(fasta)

def get_seq(site, window):
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
    # read editing sites table
    editing_sites_table=pd.read_csv(script_dir+"human_editing_sites.csv")
    # get sequece of 400bp in each side of every editing site
    seq_list = [get_seq(site,400) for i,site in editing_sites_table.iterrows()]
    # create Near_site_folding_calculator
    nc = RNAstracture.Near_site_folding_calculator("/private/common/Software/RNAstructure/RNAstructure/data_tables/","/private/apps/bin/Fold", "/private6/Projects/Yeast_Itamar_10_2022/Fold_energy/bpRNA/bpRNA.pl")
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
    print(editing_sites_table)

