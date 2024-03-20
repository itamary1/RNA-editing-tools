# SETUP ===========================================================================================
# *************************************************************************************************


# parser ------------------------------------------------------------------
library(argparse)

# create parser object
parser <- ArgumentParser()

# specify our desired options 
# by default ArgumentParser will add an help option 
parser$add_argument("-i","--input_dir", action="store", dest = "input_dir", help="Path to directory with Salmon output")
parser$add_argument("-o","--ouput", action="store", dest = "outdir", default = getwd(), help="Output directory path [.]")
#parser$add_argument("-g","--groups", action="store", dest = "groups",  help="String with names of groups, delimited by comma. Must be in the same order as input files' order", default=NULL)
parser$add_argument("-s","--salmon_suffix", action="store", dest = "salmon_suffix", default = ".quant.genes.sf", help="Unique suffix of known output files [.quant.genes.sf]")
parser$add_argument("-p","--out_suffix", action="store", dest = "out_suffix", default = "", help="Output file suffix ['']")
parser$add_argument("--wide_genes",action="store_true",dest = "output_wide_genes",default=FALSE, help="Ouput in wide format, where genes are column names")
parser$add_argument("--wide_samples",action="store_true",dest = "output_wide_samples",default=FALSE, help="Ouput in wide format, where samples are column names")
parser$add_argument("--adar",action="store",dest = "adar_pattern",default="ADAR", help="ADAR pattern to output also ADAR data only, case insensitive (matches human & mouse) [ADAR]")
parser$add_argument("--gsea", action="store_true", default = FALSE, dest = "getgsea",  help="Output normalized counts in GCT format along with phenotype file for use in GSEA analysis  [default %(default)s]")
parser$add_argument("--old_format",action="store_true",dest = "old_format",default=FALSE, help="If the files are in full old hierarchy (quant.sh within folder by dirname) than use this")
#parser$print_help()
user_args <- parser$parse_args()

# ### FOR DEBUG
# user_args$input_dir = "/private9/Projects/ImmuneCheckpointInhibitors/GSE91061_riazDataset/Salmon_1.4.0_over50b/"
# user_args$outdir = "/home/alu/fulther/Validations/"


# process -----------------------------------------------------------------
#halt execution if any of the files are missing
stopifnot(!is.null(user_args$outdir), !is.null(user_args$input_dir)) #changed to match input options

#create output directory if not exist
if (!dir.exists(user_args$outdir))
  dir.create(user_args$outdir,recursive = T)

# #get the prefix for the data
# if (is.null(user_args$prefix))
#   user_args$prefix=paste("known",gsub(",[ ]*","_", user_args$groups),"diff_0.01",sep="_")


# RUN ===========================================================================================
# *************************************************************************************************

# main --------------------------------------------------------------------
# get known files - get full file list and search for known suffix
salmon_files <- list.files(user_args$input_dir, recursive=TRUE, full.names = T)
salmon_files <- salmon_files[grep(pattern = user_args$salmon_suffix, fixed = T, x = salmon_files)]

library(data.table)
library(dplyr)
library(purrr)
library(tidyr)
library(tibble)

if (user_args$old_format) {
  salmon_out = salmon_files %>%
    set_names() %>%
    map_dfr(fread, .id = "Sample") %>%
    select(Sample, Name, TPM) %>%
    mutate(Sample = basename(dirname(Sample))) %>%
    rename(GeneSymbol=Name)
} else {
  salmon_out = salmon_files %>%
    set_names() %>%
    map_dfr(fread, .id = "Sample") %>%
    select(Sample, Name, TPM) %>%
    mutate(Sample = gsub(basename(Sample), pattern = user_args$salmon_suffix, replacement = "")) %>%
    rename(GeneSymbol=Name)
}


if (user_args$output_wide_genes) {
  # wide output genes
  salmon_out_curr <- salmon_out %>%
    pivot_wider(names_from = GeneSymbol, values_from = TPM)
  write.csv(x = salmon_out_curr, file = file.path(user_args$outdir, paste0("SalmonTPM", user_args$out_suffix, ".wideGenes.csv")), quote = F, row.names = F)
  # adar only
  salmon_out_curr %>%
    select(Sample, contains(user_args$adar_pattern, ignore.case = T)) %>%
    write.csv(file = file.path(user_args$outdir, paste0("SalmonTPM_", user_args$adar_pattern,  user_args$out_suffix, ".wideGenes.csv")), quote = F, row.names = F)
} 

if (user_args$output_wide_samples) {
  # wide output samples
  salmon_out_curr <- salmon_out %>%
    pivot_wider(names_from = Sample, values_from = TPM) 
  write.csv(x = salmon_out_curr, file = file.path(user_args$outdir, paste0("SalmonTPM", user_args$out_suffix, ".wideSamples.csv")), quote = F, row.names = F)
  # adar only
  salmon_out_curr %>%
    filter(grepl(x = GeneSymbol, user_args$adar_pattern, ignore.case = T)) %>%
    write.csv(file = file.path(user_args$outdir, paste0("SalmonTPM_", user_args$adar_pattern, user_args$out_suffix, ".wideSamples.csv")), quote = F, row.names = F)
} 

if (!user_args$output_wide_genes & !user_args$output_wide_samples) {
  # original data
  write.csv(x = salmon_out, file = file.path(user_args$outdir, paste0("SalmonTPM", user_args$out_suffix, ".csv")), quote = F, row.names = F)
  salmon_out_curr %>%
    filter(grepl(x = GeneSymbol, user_args$adar_pattern, ignore.case = T)) %>%
    write.csv(file = file.path(user_args$outdir, paste0("SalmonTPM_", user_args$adar_pattern, user_args$out_suffix, ".csv")), quote = F, row.names = F)
}
  

if(user_args$getgsea) {
  # calculate num samples + num genes from table
  numSamples = length(unique(salmon_out$Sample))
  numGenes = length(unique(salmon_out$GeneSymbol))
  normalizedCounts <- salmon_out %>% 
    pivot_wider(names_from = Sample, values_from = TPM) %>% 
    rename("NAME" = "GeneSymbol") %>%
    mutate(Description = NA) %>%
    select(NAME, Description, everything())
  
    
  header_gsea = "# 1_2"
  # write header and first line
  writeLines(paste0(header_gsea, "\n", numGenes, "\t", numSamples), 
             file.path(out,paste0(out_prefix,"normalizedCounts.gct")))
  write.table(x = normalizedCounts, file = file.path(out,paste0(out_prefix,"normalizedCounts.gct")), 
              sep = "\t", quote = F, row.names = F, append = T)
}
# ***************************************************************************************


