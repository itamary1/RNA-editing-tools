# SETUP ===========================================================================================
# *************************************************************************************************


# parser ------------------------------------------------------------------
library(argparse, quietly = T)

# create parser object
parser <- ArgumentParser()

# specify our desired options 
# by default ArgumentParser will add an help option 
parser$add_argument("-i", "--input_directory", dest="input_directory", action="store",
                    required=TRUE, help="Input directory, root")
parser$add_argument("-is", "--input_file_suffix", dest="input_file_suffix", action="store",required=TRUE, help="Suffix of files to split")
parser$add_argument("-o", "--output_dir", dest="output_dir", action="store", required=TRUE, help="Output directory of split files, will be created if does not exist")
parser$add_argument("-p", "--output_prefix", dest="output_prefix", action="store",
                    required=TRUE, help="Prefix name of output split files names")
parser$add_argument("-os", "--output_suffix", dest="output_suffix", action="store",
                    nargs="?",required=TRUE, help="Suffix of output split files names")
# parser$add_argument("-c", "--chr_list", dest="chr_list", action="store", metavar="chr_list", required=TRUE, help="List of chromosomes, speeds up the analysis, but restricts to known chromosomes")
parser$add_argument("-f", "--field", dest="field", action="store", required=TRUE, help="Field name in which to search for chromosome")
parser$add_argument("-d", "--delimiter", dest="delimiter", action="store", default=",",required=FALSE, help="Field delimiter [default %(default)s]")
parser$add_argument("--wanted_columns", dest="wanted_columns", action="store", required=FALSE,nargs="+",
                    help="List of wanted columns in output, comma separated. If none is given, all columns will be outputted")
parser$add_argument("--no_header", dest="no_header", action="store_true", default=FALSE, help="Use for headerless file, in that case columns are named R-style: V1, V2, V3, ...")
user_args <- parser$parse_args()

print(paste("Wanted columns:", user_args$wanted_columns))

# ### FOR DEBUG
# user_args <- parser$parse_args(c("-i ",
#                                 "-is StrandDerivingCountsPerRegion.csv",
#                                 "-o /home/alu/fulther/Validations",
#                                 "--output_suffix preprocessed.csv",
#                                 "--output_prefix test",
#                                 "-f GenomicRegion",
#                                 "--wanted_columns Sample AluElement Length SenseGeneCommonName SenseGeneRefSeqID AntisenseGenomicPosition SenseGenomicPosition AntisenseGeneRefSeqID AntisenseGeneCommonName A2GEditingIndex C2TEditingIndex TotalCoverageAtAllPositions MeanCoveragePerRegion NumOfA2GMismatchesNonStranded NumOfANonStranded NumOfC2TMismatchesNonStranded NumOfTNonStranded"))

# user_args$input_directory = "/private9/Projects/dsRNAProject/RNAEditingIndexOnHERegion/A2G/Plus/KidneyCortex"
# user_args$input_file_suffix = "StrandDerivingCountsPerRegion"
# user_args$output_dir = "/home/alu/fulther/Validations/"
# user_args$output_suffix = "preprocessed.csv"
# user_args$output_prefix = "test"
# user_args$field = "GenomicRegion"
# user_args$wanted_columns = c("Sample", "SenseGeneCommonName", "SenseGeneRefSeqID", "AntisenseGenomicPosition", "SenseGenomicPosition", "AntisenseGeneRefSeqID", "AntisenseGeneCommonName", "A2GEditingIndex", "C2TEditingIndex", "TotalCoverageAtAllPositions", "MeanCoveragePerRegion", "NumOfA2GMismatchesNonStranded", "NumOfANonStranded", "NumOfC2TMismatchesNonStranded", "NumOfTNonStranded")

# logging -----------------------------------------------------------------
library(log4r, quietly = T)
my_layout <- function(level, ...) {
  paste0("[",format(Sys.time()), "] ", level, "\t", ..., "\n", collapse = "")
}

logger <- logger(threshold = "DEBUG", appenders = file_appender(file.path(user_args$input_directory, format(Sys.time(), "SplitRegionsByChromosomes.%Y-%m-%dT%T.log")),
                                                                layout = my_layout))

debug(logger = logger, paste("User arguments:", user_args))



# sourcing ----------------------------------------------------------------
processPerRegionPerSampleOutput <- function(curr_file) {
  library(dplyr, quietly = T)
    curr_file <- curr_file %>%
      mutate(
        # create a "Region column
        AluElement = paste0(GenomicRegion, ":", Start, "-", End),
        Length = End - Start,
        # get complete coverage from all nucleotides
        TotalCoverageAtAllPositions = TotalCoverageAtTPositions + TotalCoverageAtCPositions + TotalCoverageAtGPositions + TotalCoverageAtAPositions) %>%
      select(Sample, AluElement, Length, ends_with("EditingIndex"), contains("Gen"),  everything(), GenomicRegion, Start, End)
  
  return(curr_file)
}


# Find files --------------------------------------------------------------
# get file list then filter by pattern
info(logger, "Finding files...")
file_list <- list.files(user_args$input_directory, recursive=TRUE, 
                        full.names = T, pattern = paste0("*", user_args$input_file_suffix))

debug(logger, paste("Found:", length(file_list), "files"))
debug(logger, paste("Files found:", paste(file_list, collapse = ", ")))


# library(data.table)
# library(dplyr)
# file_list[1] %>% 
#   fread(header=T, check.names = F, stringsAsFactors=FALSE) %>%
#   processPerRegionPerSampleOutput %>%
#   select(.data[[user_args$field]], all_of(user_args$wanted_columns)) %>%
#   group_by(.data[[user_args$field]]) %>%
#   head

#create output directory if not exist
if (!dir.exists(user_args$output_dir))
  dir.create(user_args$output_dir,recursive = T)

library(data.table, quietly = T)
library(dplyr, quietly = T, warn.conflicts = F)
library(purrr, quietly = T, warn.conflicts = F)

warn(logger, "Appending to files in wanted output directory")


# file_list %>%
#   walk(.f = compose(function (f) group_walk(.data = f, .keep = T,
#                                        .f = ~ readr::write_csv(.x, file = file.path(user_args$output_dir,
#                                                                           paste0(user_args$output_prefix, ".",
#                                                                                  unique(.x[,user_args$field]), ".",
#                                                                                  user_args$output_suffix)),
#                                                      append = T)),
#                function (f) group_by(f, .data[[user_args$field]]),
#                function (f) select(f, .data[[user_args$field]], all_of(user_args$wanted_columns)),
#                function (f) processPerRegionPerSampleOutput(f),
#                function (f) fread(f, header=T, check.names = F, stringsAsFactors=FALSE)))


for (f in file_list) {
  debug(logger, paste("Splitting", f))
  temp_file <- fread(f, header=T, check.names = F, stringsAsFactors=FALSE)
  temp_file <- processPerRegionPerSampleOutput(temp_file)
  gc()
  temp_file <- select(temp_file, .data[[user_args$field]], all_of(user_args$wanted_columns))
  gc()
  temp_file %>%
    group_by(.data[[user_args$field]]) %>%
    group_walk(.keep = T,
               .f = ~ readr::write_csv(.x, file = file.path(user_args$output_dir,
                                                            paste0(user_args$output_prefix, ".",
                                                                   unique(.x[,user_args$field]), ".",
                                                                   user_args$output_suffix)),
                                       append = T))
  gc()
}
info(logger, "Done splitting files by chromosomes")

# **********************************************************************************************************


