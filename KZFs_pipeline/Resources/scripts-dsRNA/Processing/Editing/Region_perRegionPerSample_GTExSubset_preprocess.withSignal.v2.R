###DESCRIPTION: PROCESS AEI PER REGION PER SAMPLE OUTPUT
###             SEE IF THERE IS A CERTAIN ELEMENT THAT CORRELATES TO FULL INDEX


rm(list = ls(all = TRUE))
# SETUP ===========================================================================================
# *************************************************************************************************

library(argparse, quietly = T)

# create parser object
parser <- ArgumentParser()

# specify our desired options 
# by default ArgumentParser will add an help option 
parser$add_argument("-i", "--input_directory", dest="input_directory", action="store",
                    required=TRUE, help="Input directory, root")
parser$add_argument("-o", "--output_directory", dest="output_dir", action="store", required=TRUE, help="Output directory of split files, will be created if does not exist")
parser$add_argument('-l', '--log_path', dest='log_path', action='store', default="", help='Log file, default is to create in input dir')
parser$add_argument("-c", "--columns", dest="columns", action="store", required=TRUE, nargs="+", help="Column names for file")
parser$add_argument("-g", "--group_file", dest="group_file", action="store", required=TRUE, help="CSV with sample to group. Will pool samples according to given groups, sample count per group will be calculated according to this file. Format: Sample,Group")
parser$add_argument("-v", "--group_vars", dest="grouping_variables", action="store", required=TRUE, nargs="+", help="Variables to group by")
parser$add_argument("-p", "--prefix", dest="prefix", action="store", default="StrandDerivingCountsPerRegion", help="Prefix of files to group (before .chr?.)")
parser$add_argument("-s", "--suffix", dest="suffix", action="store", default="preprocessed.csv", help="Suffix of files to group (after .chr?.)")
parser$add_argument("-op", "--output_prefix", dest="output_prefix", action="store", default="EditingPerRegionPerSample", nargs="+", help="Prefix of output files (before .chr?.)")
parser$add_argument("-os", "--output_suffix", dest="output_suffix", action="store", default="pooled.withSignal.csv", nargs="+", help="Suffix of output files (after .chr?.)")
user_args <- parser$parse_args()

# logging -----------------------------------------------------------------
library(log4r, quietly = T, warn.conflicts = F)
library(data.table, quietly = T)
library(dtplyr, quietly = T)
library(dplyr, warn.conflicts = FALSE, quietly = T)

my_layout <- function(level, ...) {
  paste0("[",format(Sys.time()), "] ", level, "\t", ..., "\n", collapse = "")
}

logger <- logger(threshold = "DEBUG", appenders = file_appender(file = if_else(user_args$log_path == "", 
                                                          true = file.path(user_args$input_directory, paste0("Pool", user_args$prefix,format(Sys.time(), "%Y-%m-%dT%T.log"))), 
                                                          false = user_args$log_path),
                                           layout = my_layout))

debug(logger = logger, paste("User arguments:", user_args))

# General -----------------------------------------------------------------

### Per region per sample, Alu index (via Hillel's pipeline file, summary files: "EditingIndex.csv")
EI_perRegionPerSample_files = list.files(user_args$input_directory, 
                                   pattern = paste0(user_args$prefix, ".chr.*", user_args$suffix),
                                   full.names = T, recursive = F)
debug(logger = logger, paste("Files found:", paste(EI_perRegionPerSample_files, collapse = ", ")))

### Output directory
outdir = user_args$output_dir
# out_boxplots = "boxplots"
# out_scatterplots = "scatterplots"
# out_correlations = "correlations"
# out_barplots = "barplots"
# out_hist = "histograms"
out_stats = "stats"
# out_stats_plus = file.path(out_stats, "Plus")
# out_stats_minus = file.path(out_stats, "Minus")
# out_plots = "plots"
# out_bed = "bed_files"
# out_bed_lists = "lists"
# out_bed_elements = "elements"

# PREPROCESSING ===================================================================================
# *************************************************************************************************

# create directories ------------------------------------------------------
#create paths, change
# out_plots = file.path(outdir, out_plots)
out_stats = file.path(outdir, out_stats)
# out_stats_plus = file.path(outdir, out_stats_plus)
# out_stats_minus = file.path(outdir, out_stats_minus)
# out_boxplots = file.path(out_plots, out_boxplots)
# out_barplots = file.path(out_plots, out_barplots)
# out_scatterplots = file.path(out_plots, out_scatterplots)
# out_correlations = file.path(out_plots, out_correlations)
# out_hist = file.path(out_plots, out_hist)
# out_bed = file.path(outdir, out_bed)
# out_bed_lists = file.path(out_bed, out_bed_lists)
# out_bed_elements = file.path(out_bed, out_bed_elements)
#create directories, if do not exist
# dir.create(out_boxplots, recursive = T)
# dir.create(out_barplots, recursive = T)
# dir.create(out_scatterplots, recursive = T)
# dir.create(out_correlations, recursive = T)
# dir.create(out_hist, recursive = T)
dir.create(out_stats, recursive = T)
# dir.create(out_bed_lists, recursive = T)
# dir.create(out_bed_elements, recursive = T)


# load & process SRA -------------------------------------------------------------
info(logger, "Reading classification file and calculating group size")
sra = read.csv(file=user_args$group_file, header=T, check.names = F, stringsAsFactors=FALSE,
               col.names = c("Sample", "Group")) %>%
  # add group count
  group_by(Group) %>% 
  add_tally(name = "TotalSamplesPerGroup")  

# colors -------------------------------------------------------------------
# GTEx_colors <- read.delim(GTEx_colors, header = T, check.names = F, stringsAsFactors = F) 


# grouping ----------------------------------------------------------------
column_names=user_args$columns
group_by_vars=user_args$grouping_variables

# ANALYSIS =======================================================================================
# *************************************************************************************************
info(logger, "Pooling each chromosome")

# pooling of core date ---------------------------------------------------------------
for (f in c(EI_perRegionPerSample_files)) {
  # get chromosome number
  chromosome = gsub(basename(f), pattern = paste0(user_args$prefix, ".|.", user_args$suffix), replacement = "")
  debug(logger, paste("Loading file", f))
  EI_perRegionPerSample <- fread(file = f, header = F, stringsAsFactors = F,
                                 col.names =  column_names) %>%
    inner_join(sra) %>%
    as.data.table()
  
  # change suffix
  output_suffix = paste0(".", chromosome, ".")
  debug(logger, paste("Beginning pooling analysis for", chromosome))
  
  ## pool -------------------------------------------
  # pool by element
  EI_perRegionPerSample_pooled = EI_perRegionPerSample %>%
    # lazy_dt() %>%
    group_by(across(all_of(group_by_vars))) %>%
    summarise(across(c(Length:AntisenseGeneCommonName, TotalSamplesPerGroup), first),
              across(matches("IndexedCanonicalOf[ATCG]2[ATCG]|IndexedMismatchesOf[ATCG]2[ATCG]|TotalCoverageAtAllPositions"), sum),
              across(matches("NumOfIndexed(Mismatches|Overall)SitesOf[ACGT]2[ACGT]"), list(Min=min, Max=max), .names = "{.fn}{.col}"),
              SampleCount = n()) %>%
    mutate(MeanCoveragePerRegion = TotalCoverageAtAllPositions / Length,
           
           A2GEditingIndex = 100 * IndexedMismatchesOfA2G /  (IndexedMismatchesOfA2G + IndexedCanonicalOfA2G),
           A2CEditingIndex = 100 * IndexedMismatchesOfA2C /  (IndexedMismatchesOfA2C + IndexedCanonicalOfA2C),
           A2TEditingIndex = 100 * IndexedMismatchesOfA2T /  (IndexedMismatchesOfA2T + IndexedCanonicalOfA2T),
           
           C2AEditingIndex = 100 * IndexedMismatchesOfC2A /  (IndexedMismatchesOfC2A + IndexedCanonicalOfC2A),
           C2GEditingIndex = 100 * IndexedMismatchesOfC2G /  (IndexedMismatchesOfC2G + IndexedCanonicalOfC2G),
           C2TEditingIndex = 100 * IndexedMismatchesOfC2T /  (IndexedMismatchesOfC2T + IndexedCanonicalOfC2T),
           
           G2AEditingIndex = 100 * IndexedMismatchesOfG2A /  (IndexedMismatchesOfG2A + IndexedCanonicalOfG2A),
           G2CEditingIndex = 100 * IndexedMismatchesOfG2C /  (IndexedMismatchesOfG2C + IndexedCanonicalOfG2C),
           G2TEditingIndex = 100 * IndexedMismatchesOfG2T /  (IndexedMismatchesOfG2T + IndexedCanonicalOfG2T),
           
           T2AEditingIndex = 100 * IndexedMismatchesOfT2A /  (IndexedMismatchesOfT2A + IndexedCanonicalOfT2A),
           T2CEditingIndex = 100 * IndexedMismatchesOfT2C /  (IndexedMismatchesOfT2C + IndexedCanonicalOfT2C),
           T2GEditingIndex = 100 * IndexedMismatchesOfT2G /  (IndexedMismatchesOfT2G + IndexedCanonicalOfT2G),
           
           AvgMeanCoveragePerRegionInGroup = MeanCoveragePerRegion / TotalSamplesPerGroup) %>%
    select(AluElement, Length, contains("EditingIndex"), AvgMeanCoveragePerRegionInGroup, MeanCoveragePerRegion, 
           SampleCount,
           SenseGeneCommonName, AntisenseGeneCommonName, everything()) %>%
    as.data.table()
  outfile = file.path(out_stats, paste0(user_args$output_prefix, output_suffix, user_args$output_suffix))
  debug(logger, paste("Writing", f, "into", outfile))
  fwrite(EI_perRegionPerSample_pooled, outfile, quote = F, row.names = F, scipen = 999)
  
  debug(logger, paste("Summing sample counts for", chromosome))
  ## sample and Group information -------------------------------------------
  # pool by element
  EI_perRegionPerSample_pooled_count = EI_perRegionPerSample %>%
    # lazy_dt() %>%
    group_by(across(all_of(group_by_vars))) %>%
    # count number of samples
    summarise(SampleCount = n(),
              # at least A2G MM = edited sample
              SampleCountEdited = sum(IndexedMismatchesOfA2G > 0, na.rm = T)) %>%
    mutate(EditStatus = if_else(SampleCountEdited > 0, "Edited", "NotEdited")) %>%
    group_by(AluElement, EditStatus) %>%
    summarise(GroupCount = n(),
              # GroupCountEdited = sum(SampleCountEdited > 0, na.rm = T),
              Groups = paste0(unique(Group), collapse = ";"),
              across(SampleCount:SampleCountEdited, sum)) %>%
    select(-SampleCountEdited) %>%
    tidyr::pivot_wider(id_cols = AluElement,
                       names_from = EditStatus, 
                       values_from = c(contains("Count"), contains("Group"))) %>%
    tidyr::replace_na(replace = list(GroupCount_Edited = as.integer(0), SampleCount_Edited = as.integer(0),
                                     GroupCount_NotEdited = as.integer(0), SampleCount_NotEdited = as.integer(0))) %>%
    as.data.table()
  outfile = file.path(out_stats, paste0("CountDetails", user_args$output_prefix, output_suffix, user_args$output_suffix))
  debug(logger, paste("Writing", f, "into", outfile))
  fwrite(EI_perRegionPerSample_pooled_count, outfile, quote = F, row.names = F, scipen = 999)
}


info(logger, "Done pooling all chromosomes")


# **********************************************************************************************************


