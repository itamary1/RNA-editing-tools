#!/bin/sh
script_dir=$(dirname ${0})
genome=${script_dir}"/first10000_sacCer3.fa"
editing_sites_table=${script_dir}"/editing_sites.csv"
data_tables="RNAstructure/data_tables/"
FOLD="Fold"
bpRNA="bpRNA.pl"
out_table_path=${script_dir}"/editing_sites_with_substructurs.csv"

python Substructures_Calculator.py $genome $editing_sites_table $data_tables $FOLD $bpRNA $out_table_path

