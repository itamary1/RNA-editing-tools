library(dplyr);library(tidyr) ;library(data.table);library(argparse)

# create parser object
parser <- ArgumentParser(description='for creating rid file from bed6 containing gene_name at column 4. rid file is region-id csv file')

parser$add_argument("-i", action="store", dest = "infile", help="Input file, bed6 format containing gene_name at column 4")
parser$add_argument("-o", action="store", dest = "outfile",  help="Output RID file")

user_args <- parser$parse_args()

fread(user_args$infile, header = F, stringsAsFactors = F, col.names = c("chr", "start", "end", "gene", "score", "strand")) %>% separate_rows(gene, sep = ",") %>% mutate(region = paste0(chr, ":", start,  "-", end)) %>% select(region,gene) %>% fwrite(user_args$outfile, quote = F, row.names = F, scipen = 999)