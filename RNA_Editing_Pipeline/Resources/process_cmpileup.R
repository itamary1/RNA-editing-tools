# SETUP ===========================================================================================
# *************************************************************************************************
rm(list = ls(all = TRUE))
library(argparse, quietly = T, warn.conflicts = F)

# create parser object
parser <- ArgumentParser()

# specify our desired options 
# by default ArgumentParser will add an help option 
parser$add_argument("-i", "--input_dir", dest="input_dir", action="store", nargs='+',
                    required=TRUE, help="Input directory with cmpileup files")
parser$add_argument('-o', '--output_dir',dest="output_dir", default="", help='Path of processed results directory, default is to create in input dir', required=FALSE)
parser$add_argument('-l', '--log_path', dest='log_path', action='store', metavar='log_path', default="", help='Log file, default is to create in input dir')
parser$add_argument('-s', '--suff_cmpileup', dest="suff_cmpileup", default="_mpileup.cmpileup", help='cmpileup files suffix')
parser$add_argument('-c', '--complement', dest="complement", action="store_true", default=FALSE, help='Create also a complemented version')
parser$add_argument("-p", "--process_num", dest="process_num", action="store", type="integer", default=4, help="Number of processes to run with")
parser$add_argument("-f", "--force", dest="force", action="store_true", default=FALSE, help="Force overwrite")
user_args <- parser$parse_args()

### Output directory
indir = user_args$input_dir
outdir = user_args$output_dir
force_run = user_args$force
num_processes = user_args$process_num
cmpileup_suffix = user_args$suff_cmpileup
complement_files = user_args$complement
if(user_args$log_path == "") {
  log_dir = indir
} else {
  log_dir = user_args$log_path
}

# logging -----------------------------------------------------------------
library(log4r, quietly = T, warn.conflicts = F)
library(data.table, quietly = T, warn.conflicts = F)
library(dplyr, quietly = T, warn.conflicts = F)

my_layout <- function(level, ...) {
  paste0("[",format(Sys.time()), "]\t", level, "\t", Sys.getpid(), "\t", ..., "\n", collapse = "")
}

logger <- logger(threshold = "DEBUG", appenders = file_appender(file.path(log_dir, format(Sys.time(), "ProcessCMpileup.%Y-%m-%dT%T.log")), layout = my_layout))

debug(logger = logger, paste("User arguments:", user_args))

# Find Files -----------------------------------------------------------------
input_files_list = list.files(path = indir, pattern = paste0(cmpileup_suffix, "$"), recursive = T, full.names = T)
debug(logger, paste("Files:", input_files_list))

# helper functions --------------------------------------------------------
# calculate index for single combination
index <- function(df, mismatch, canonical) {
  canonical_base = df %>% ungroup %>% select({{canonical}}) %>% colnames()
  mutate(df, "{{canonical}}2{{mismatch}}EditingIndex" := if_else(Reference == canonical_base, 
                                                                 true = 100 * {{mismatch}} / ({{mismatch}} + {{canonical}}),
                                                                 false = NA_real_))
}

# calculate all possible indexes
process <- function(df) {
  mutate(.data = df, 
         AGCoveragePerSite = A + G) %>%
    index(canonical = A, mismatch = C) %>%
    index(canonical = A, mismatch = G) %>%
    index(canonical = A, mismatch = `T`) %>%
    index(canonical = C, mismatch = A) %>%
    index(canonical = C, mismatch = G) %>%
    index(canonical = C, mismatch = `T`) %>%
    index(canonical = G, mismatch = A) %>%
    index(canonical = G, mismatch = C) %>%
    index(canonical = G, mismatch = `T`) %>%
    index(canonical = `T`, mismatch = A) %>%
    index(canonical = `T`, mismatch = C) %>%
    index(canonical = `T`, mismatch = G) %>%
    mutate(across(ends_with("EditingIndex"), ~na_if(., NA)))
}

complement <- function(df) {
  df = df %>%
    # first swap count columns names to complementary nt
    rename(A = `T`,C = G,G = C, `T` = A) %>%
    # swap reference to complementary nt - will now represent the *indexed* canonical base
    mutate(Reference = case_when(Reference == "A" ~ "T",
                                 Reference == "C" ~ "G",
                                 Reference == "G" ~ "C",
                                 Reference == "T" ~ "A",
                                 TRUE ~ stringr::str_to_upper(Reference)))
  return(df)
}

format_columns <- function(df) {
  # match processed cmpileup format
  mutate(.data = df,
         Region = paste0(RegionChr, ":", RegionStart, "-", RegionEnd),
         PositionChr = RegionChr,
         PositionStart = PositionInRegion - 1,
         PositionEnd = PositionInRegion,
         Reference = stringr::str_to_upper(Reference)) %>%
    select(PositionChr, PositionStart, PositionEnd,
           Region,
           Reference,
           TotalCoverage, A, C, G, `T`, UnrecognizedBases, LowQualityBases)
}

# PROCESS =========================================================================================
# *************************************************************************************************


# process -------------------------------------------------------------
info(logger, "Processing file")

process_and_complement_files <- function(f) {
  ## plus -----
  # process files as is - assume all regions are on + strand and A2G is the wanted MM
  processed_filename=paste0(f, ".processed.csv")
  # if exists and should not force overwrite - skip
  if (file.exists(processed_filename) & !force_run) {
    warn(logger, paste("File already processed on plus strand, delete or force: ", processed_filename))
  } else {
    debug(logger, paste("Processing", f))
    fread(f, stringsAsFactors = F, header = F, col.names = c("RegionChr", "RegionStart", "RegionEnd", 
                                                             "PositionInRegion", "Reference",
                                                             "TotalCoverage", 
                                                             "A", "C", "G", "T", 
                                                             "UnrecognizedBases", "LowQualityBases")) %>% 
      # process reference and column names 
      format_columns %>%
      # add index
      process %>%
      fwrite(processed_filename, quote = F, row.names = F, scipen = 999)
    # inform to ensure writing is finished
    debug(logger, paste("Finished writing", processed_filename))
    gc()
  }
  
  if (complement_files)  {
    # minus -----
    # swap files to match complementary strand
    # then process new files - assume all regions are on - strand and T2C is the wanted MM
    processed_filename=paste0(f, ".processed.complemented.csv")
    if (file.exists(processed_filename) & !force_run) {
      warn(logger, paste("File already processed on minus strand, delete or force: ", processed_filename))
    } else {
      debug(logger, paste("Processing and complementing", f))
      fread(f, stringsAsFactors = F, header = F, col.names = c("RegionChr", "RegionStart", "RegionEnd", 
                                                               "PositionInRegion", "Reference",
                                                               "TotalCoverage", 
                                                               "A", "C", "G", "T", 
                                                               "UnrecognizedBases", "LowQualityBases")) %>% 
        # process reference and column names 
        format_columns %>%
        # complement - take T reference and swap bases
        complement %>%
        # add index
        # process as plus after complementation
        process %>%
        fwrite(processed_filename, quote = F, row.names = F, scipen = 999)
      # inform to ensure writing is finished
      debug(logger, paste("Finished writing", processed_filename))
      gc()
    } 
  }
}


# main  -------------------------------------------------------------------
# set processing plan to multicore and set number of cores
options(future.globals.maxSize= 70000 * 1024 ^ 2)
future::plan(future::multicore, workers = num_processes)


# call processing
furrr::future_walk(.x = input_files_list,
                   .f = ~ process_and_complement_files(f = .x))


info(logger, "Done processing all CMpileup files")


# **************************************************************************

