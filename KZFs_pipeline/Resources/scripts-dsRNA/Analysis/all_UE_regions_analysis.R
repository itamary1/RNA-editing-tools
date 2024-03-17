
rm(list = ls(all = TRUE))
# SETUP ===========================================================================================
# *************************************************************************************************

library(argparse, quietly = T)

# create parser object
parser <- ArgumentParser()

# specify our desired options 
# by default ArgumentParser will add an help option 
parser$add_argument("-r", "--pooled_regions", dest="pooled_region_editing", action="store",
                    required=TRUE, help="File with pooled editing as created by the `regions` program")
parser$add_argument("-s", "--pooled_sites", dest="pooled_site_editing", action="store", required=TRUE, help="File with pooled editing per region from site information as created by the `sites` program")
parser$add_argument("-c", "--editing_count", dest="editing_count", action="store", required=TRUE, nargs="+", help="File with edited and total sample count per region as created by `regions` program")
parser$add_argument("-g", "--group_file", dest="group_file", action="store", required=TRUE, help="CSV with sample to group. Will pool samples according to given groups, sample count per group will be calculated according to this file. Format: Sample,Group")
parser$add_argument("-o", "--output_directory", dest="output_dir", action="store", required=TRUE, help="Output directory of split files, will be created if does not exist")
parser$add_argument('-l', '--log_path', dest='log_path', action='store', default="", help='Log file, default is to create in input dir')
#parser$add_argument("-v", "--group_vars", dest="grouping_variables", action="store", required=TRUE, nargs="+", help="Variables to group by")
#parser$add_argument("-p", "--prefix", dest="prefix", action="store", default="StrandDerivingCountsPerRegion", help="Prefix of files to group (before .chr?.)")
#parser$add_argument("-op", "--output_prefix", dest="output_prefix", action="store", default="EditingPerRegionPerSample", nargs="+", help="Prefix of output files (before .chr?.)")
#parser$add_argument("-os", "--output_suffix", dest="output_suffix", action="store", default="pooled.withSignal.csv", nargs="+", help="Suffix of output files (after .chr?.)")
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


# External Scripts --------------------------------------------------------

# General -----------------------------------------------------------------

### SRA table
sra = user_args$group_file

# sample count
sample_count = sra %>%
  # add group count
  group_by(Group) %>% 
  add_tally(name = "TotalSamplesPerGroup")  

# preprocessed AEI files
#"/private9/Projects/dsRNAProject/Summary_AluEditingIndexWithHEReads/stats/AluRegionsPerRegionPerSample.pooled.withSignal.csv"
Alu_EI_perRegionPerSample_pooled = user_args$pooled_region_editing

# count AEI information files
#"/private9/Projects/dsRNAProject/Summary_AluEditingIndexWithHEReads/stats/CountDetailsAluRegionsPerRegionPerSample.pooled.withSignal.csv"
count_details = user_args$editing_count

# pooled site information
#"/private9/Projects/dsRNAProject/Summary/AluEditingIndexWithHEReads/AluEditedRegionsSitesInformation.pooled.csv"
Alu_nt_count_per_ES_pooled_aggregated = user_args$pooled_site_editing

### Params
# total number of samples
num_subset_samples = sum(fread(sample_count)$TotalSamplesPerTissueAndCohort)
num_tissues = length(unique(fread(sample_count)$Tissue))
num_tissues_cohort = length(unique(paste0(fread(sample_count)$Tissue, fread(sample_count)$Cohort)))


### Output directory
outdir = user_args$output_dir

# PREPROCESSING ===================================================================================
# *************************************************************************************************

# create directories ------------------------------------------------------
out_plots = file.path(outdir)


# load & process SRA -------------------------------------------------------------
sra = read.csv(file=sample_count,header=T, check.names = F, stringsAsFactors=FALSE)


# sample count ------------------------------------------------------------
sample_count = sra %>%
  # add group count
  group_by(Group) %>% 
  add_tally(name = "TotalSamplesPerGroup")  


# Count information ----------------------------------------------------------
debug(logger, "Read and process region count details")
count_details = count_details %>% 
  purrr::map_dfr(fread, header = T, stringsAsFactors = F, 
                 sep =",", sep2 = ";", fill= T) %>%
  rename(!!region_col := .data[[original_region_col]])


# load and process ES editing index ---------------------------------------
debug(logger, "Read and process index per site")
Alu_nt_count_per_ES_pooled_aggregated = fread(file=Alu_nt_count_per_ES_pooled_aggregated,header=T, 
                                              check.names = F, stringsAsFactors=FALSE) %>%
  select(.data[[region_col]], contains("EditingSitesInRegion"), ends_with("EditingIndexNoSNP"), 
         starts_with("NumNonSNPASitesInRegion"), starts_with("%NonSNPASitesInRegion"), 
         TotalIndexedSitesInRegion_ReferenceBaseA, starts_with("TotalSNPSitesInRegion"))


# load and process AEI -------------------------------------------------------------
editing_levels_factors =  c("[0-5)","[5-10)","[10-15)","[15-20)", "[20-25)","25+")
coverage_range_factors = c("[0-0.25)", "[0.25-0.5)", "[0.5-1)", "[1-10)","10+")

#**CHANGED** or #**ADDED**
#*Added the filtering to create distinct groups of coverage
#*Changed the `MinAvgMeanCov` to be a [) format instead of minimal coverage
#*Loading then filtering
timestamp(prefix = "", suffix = "\t DEBUG Read and process RNA editing index")
Alu_EI_perRegionPerSample_pooled = fread(file = Alu_EI_perRegionPerSample_pooled, header = T, stringsAsFactors = F) %>%
  rename(!!region_col := .data[[original_region_col]]) %>%
  # genomic annotations are available from index itself
  # add repeat annotation information
  inner_join(annotation_repeats %>%
               select(-OverlapFrac)) %>%
  # add count information
  inner_join(count_details, by = region_col) %>%
  mutate(across(contains("SampleCount"), ~ .x *100 / TotalSamples, .names = "{.col}%"),
         across(contains("TissueCount"), ~ .x *100 / 47, .names = "{.col}%")) %>%
  # add editing sites information
  inner_join(Alu_nt_count_per_ES_pooled_aggregated) %>%
  mutate(A2GEditingIndexCalcDiff = abs(A2GEditingIndex - A2GEditingIndexNoSNP)) %>%
  inner_join(Alu_EI_perRegionPerSample_pooled_otherMM) %>%
  # act on row
  rowwise() %>%
  # sum all non-A2G editing index, for full and non-SNP version separately, DO NOT include NA
  # must be done here as files do not include C2T column (already calculated in major calculation)
  mutate(NonA2GEditingIndexSum = sum(c_across(matches("[CGT]2[ACGT]EditingIndex$") | matches("A2[CT]EditingIndex$")), 
                                     na.rm = T),
         NonA2GEditingIndexSumNoSNP = sum(c_across(matches("[CGT]2[ACGT]EditingIndexNoSNP$") | matches("A2[CT]EditingIndexNoSNP$")), 
                                          na.rm = T),
         # now fix index sum for overlapping regions
         NonA2GEditingIndexSumOverlapConsidered = if_else(Overlap == "Overlapping", 
                                                          true = NonA2GEditingIndexSum - T2CEditingIndex,
                                                          false = NonA2GEditingIndexSum),
         NonA2GEditingIndexSumNoSNPOverlapConsidered = if_else(Overlap == "Overlapping", 
                                                               true = NonA2GEditingIndexSumNoSNP - T2CEditingIndexNoSNP,
                                                               false = NonA2GEditingIndexSumNoSNP)) 

# iteratively classifies regions according to the following:
# For each x=5,10,15, ...
# Include regions with index >=x & at least 3 non-SNP sites with editing>=x and >=2G, unless they belong to a higher x
get_iterative_editing_classes <- function(df, region_col, global_editing_col, site_editing_col, site_cov_col, 
                                          ordered_min_editing_vec, min_site_num, min_G_cov, wanted_col) {
  # initialize
  classified_regions = c()
  top_limit = NA
  # for each limit
  for (min_editing in rev(ordered_min_editing_vec)) {
    # get current column to work with
    curr_filter_col = paste0("NumNonSNPASitesInRegionWith_", site_editing_col, "AtLeast", min_editing, "_Num", site_cov_col, "AtLeast", min_G_cov)
    
    # filter regions that pass current criteria and have not yet been classified
    if (is.na(top_limit)) {
      new_df = df %>%
        # annotate
        mutate(!!wanted_col := if_else(.data[[global_editing_col]] >= min_editing & 
                                         .data[[curr_filter_col]] >= min_site_num & 
                                         !.data[[region_col]] %in% classified_regions,
                                       true = paste0(min_editing, "+"),
                                       false = NA_character_))
    } else {
      new_df = new_df %>%
        # annotate
        mutate(!!wanted_col := if_else(.data[[global_editing_col]] >= min_editing & 
                                         .data[[curr_filter_col]] >= min_site_num & 
                                         !.data[[region_col]] %in% classified_regions,
                                       true = paste0("[", min_editing, "-", top_limit, ")"),
                                       false = .data[[wanted_col]]))
    }
    
    # update minimal editing parameter
    top_limit = min_editing
    # update list of classified regions
    classified_regions = new_df %>% filter(!is.na(.data[[wanted_col]])) %>% pull(.data[[region_col]])
  }
  
  # change NA to better value
  new_df = new_df %>%
    mutate(!!wanted_col := tidyr::replace_na(.data[[wanted_col]], "Low Editing"))
  
  return(new_df)
}


get_catagories_classes <- function(given_data, given_col, wanted_col, pairs) {
  catagories_data_classes = data.frame()
  for(pair in pairs[1:length(pairs)-1]) {
    bottom_limit = pair[1]
    top_limit = pair[2]
    curr_data <- given_data %>%
      as_tibble() %>%
      filter(.data[[given_col]] >= bottom_limit, .data[[given_col]] < top_limit) %>%
      mutate(!!wanted_col := paste0("[", bottom_limit, "-", top_limit, ")"))
    catagories_data_classes <- rbind(catagories_data_classes, curr_data)
  }
  
  bottom_limit=pairs[[length(pairs)]][1]
  curr_data <- given_data %>%
    as_tibble() %>%
    filter(.data[[given_col]] >= bottom_limit) %>%
    mutate(!!wanted_col := paste0(bottom_limit, "+"))
  catagories_data_classes <- rbind(catagories_data_classes, curr_data)
  
  return(catagories_data_classes)
}

# add editing catagories
timestamp(prefix = "", suffix = "\t DEBUG Add editing catagories")
Alu_EI_perRegionPerSample_pooled = get_catagories_classes(given_data = Alu_EI_perRegionPerSample_pooled, 
                                                          given_col = editing_col,
                                                          wanted_col = "GlobalEditingPcntClass",
                                                          pairs = list(c(0,5), c(5,10),
                                                                       c(10, 15), c(15,20),
                                                                       c(20,25), c(25,Inf))) 
Alu_EI_perRegionPerSample_pooled = get_catagories_classes(given_data = Alu_EI_perRegionPerSample_pooled, 
                                                          given_col = "A2GEditingIndexNoSNP",
                                                          wanted_col = "GlobalEditingPcntNoSNPClass",
                                                          pairs = list(c(0,5), c(5,10),
                                                                       c(10, 15), c(15,20),
                                                                       c(20,25), c(25,Inf))) 
timestamp(prefix = "", suffix = "\t DEBUG Add coverage catagories")
Alu_EI_perRegionPerSample_pooled = get_catagories_classes(given_data = Alu_EI_perRegionPerSample_pooled,
                                                          given_col  = cov_col,
                                                          wanted_col = "AverageMeanCovarageClass",
                                                          pairs = list(c(0, 0.25), c(0.25, 0.5),
                                                                       c(0.5, 1), c(1, 10),
                                                                       c(10, Inf)))
Alu_EI_perRegionPerSample_pooled = get_catagories_classes(given_data = Alu_EI_perRegionPerSample_pooled,
                                                          given_col  = "MeanCoveragePerRegion",
                                                          wanted_col = "MeanCovarageClass",
                                                          pairs = list(c(0, 100), c(100, 1000), 
                                                                       c(1000, 5000), c(5000, 10000),
                                                                       c(10000, Inf)))

timestamp(prefix = "", suffix = "\t DEBUG Add iterative new editing catagories")
Alu_EI_perRegionPerSample_pooled = get_iterative_editing_classes(df = Alu_EI_perRegionPerSample_pooled, 
                                                                 global_editing_col = "A2GEditingIndexNoSNP",
                                                                 site_editing_col = "A2GEditingIndex",
                                                                 site_cov_col = "G",
                                                                 region_col = region_col,
                                                                 min_site_num = 3, 
                                                                 min_G_cov = 2, 
                                                                 wanted_col = "EditingPcntNoSNPClass", 
                                                                 ordered_min_editing_vec = seq(5, 30, by = 5))


timestamp(prefix = "", suffix = "\t DEBUG Write output")
fwrite(Alu_EI_perRegionPerSample_pooled, file.path(outdir, paste0("AluRegionsPerRegionPerSample.AllInfo.csv")), quote = F, row.names = F, scipen = 999)
fwrite(Alu_EI_perRegionPerSample_pooled %>%
         select(.data[[region_col]], Length, Strand, Overlap, 
                matches("^A2GEditingIndex"), matches("^(Non){0,1}A2GEditingIndex"),
                AvgMeanCoveragePerRegion, MeanCoveragePerRegion,
                SampleCount, `SampleCount%`, matches("Count_Edited[%]*"), Tissues_Edited, 
                ends_with("GenomicLocation"), NumRepeats, RepeatFamily:RepeatStrandSequence, 
                matches("UEClusters.*Count"), UEClustersTissueInfo, 
                matches("^[ATCG]2[ATCG]EditingIndex"), matches("^(Non){0,1}[ATCG]2[ATCG]EditingIndex"), 
                TotalIndexedSitesInRegion_ReferenceBaseA,TotalSNPSitesInRegion_ReferenceBaseA, TotalSNPSitesInRegion, ends_with("Class"),
                starts_with("NumNonSNPASites") & ends_with("_NumGAtLeast2")), 
       file.path(outdir, paste0("AluRegionsPerRegionPerSample.AllInfo.short.csv")), quote = F, row.names = F, scipen = 999)

timestamp(prefix = "", suffix = "\t INFO Done")



# *************************************************************************************










