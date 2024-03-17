# create argsparse
ARGPARSE_DESCRIPTION="Creates the commands for a GEO dataset"      # this is optional
source /private/common/Software/BashLibs/argparse-bash/argparse.bash || exit 1
argparse "$@" <<EOF || exit 1
parser.add_argument('-i', '--input_file',type=str, help='Input CSV of summed cmpileup', required=True)
parser.add_argument('-o', '--output_file',type=str, help='Output file of non-SNP sites', required=True)
parser.add_argument('-os', '--snp_output_file',type=str, help='Output file of SNP sites', required=True)
parser.add_argument('--snps', type=str, help='Path of SNP BED file', default="/private/common/Software/AEI/RNAEditingIndex1.1/RNAEditingIndexer/Resources/SNPs/HomoSapiens/ucscHg38CommonGenomicSNPs150.bed.gz")
EOF

echo "Running with parameters: -i" $INPUT_FILE "-o" $OUTPUT_FILE "-os" $SNP_OUTPUT_FILE "--snps" $SNPS

# use A as reference, take columns of A and G
# no SNP
awk 'BEGIN{FS=","; OFS="\t"} NR > 1 {print $1, $2, $3, $4, $5, $6, $7, $8, $9, $10, $11, $12}' $INPUT_FILE | bedtools.2.30 intersect -a stdin -b $SNPS -v | awk 'BEGIN{OFS=","; print "PositionChr", "PositionStart", "PositionEnd", "Region", "Reference", "TotalCoverage", "A", "C", "G", "T", "UnrecognizedBases", "LowQualityBases"} {print $1, $2, $3, $4, $5, $6, $7, $8, $9, $10, $11, $12}'> $OUTPUT_FILE

# with SNP
awk 'BEGIN{FS=","; OFS="\t"} NR > 1 {print $1, $2, $3, $4, $5, $6, $7, $8, $9, $10, $11, $12}' $INPUT_FILE | bedtools.2.30 intersect -a stdin -b $SNPS -wa -wb | awk 'BEGIN{OFS=","; print "PositionChr", "PositionStart", "PositionEnd", "Region", "Reference", "TotalCoverage", "A", "C", "G", "T", "UnrecognizedBases", "LowQualityBases", "SNPPosition", "SNPStrand", "SNPRef", "SNPVariant"} {print $1, $2, $3, $4, $5, $6, $7, $8, $9, $10, $11, $12, $13 ":" $14 "-" $15, $16, $17, $18}' > $SNP_OUTPUT_FILE
