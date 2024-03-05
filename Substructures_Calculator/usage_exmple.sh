#!/bin/sh
script_dir=$(dirname ${0})
genome=${script_dir}"/first10000_sacCer3.fa"
editing_sites_table=${script_dir}"/editing_sites.csv"
data_tabels="/private/common/Software/RNAstructure/RNAstructure/data_tables/"
FOLD="/private/apps/bin/Fold"
bpRNA="/private6/Projects/Yeast_Itamar_10_2022/Fold_energy/bpRNA/bpRNA.pl"
out_table_path=${script_dir}"/editing_sites_with_substructurs.csv"

python /home/alu/twerski/Public/RNA-editing-tools/Substructures_Calculator/Substructures_Calculator.py $genome $editing_sites_table $data_tabels $FOLD $bpRNA $out_table_path

