# SETUP ===========================================================================================
# *************************************************************************************************
rm(list = ls(all = TRUE))
library(argparse, quietly = T)

# create parser object
parser <- ArgumentParser()

# specify our desired options 
# by default ArgumentParser will add an help option 
parser$add_argument("-ins", "--input_no_snps", dest="input_no_snps", action="store",
                    required=TRUE, help="Input file without SNPs")
parser$add_argument("-ios", "--input_only_snps", dest="input_only_snps", action="store",
                    required=TRUE, help="Input file only with SNPs")
parser$add_argument("-o", "--output_directory", dest="output_dir", action="store", required=TRUE, help="Output directory of processed files, will be created if does not exist")
parser$add_argument("-g", "--group_name", dest="group_name", action="store", required=TRUE, help="Name of group")
parser$add_argument("-n", "--group_size", dest="current_group_size", type="integer", action="store", required=TRUE, help="Number of samples in group")
parser$add_argument("-rid", "--region_id_file", dest="region_id_file", default="", action="store", help="A CSV file to replace region coordinates as grouping variable during pooling (only affects last step). Must conform to the format 'Region,NewRegionID', to ensure proper pooling. This allows calculating pooled editing over sites from non-consecutive genomic regions, by renaming them to a single region ID. Note that a left join is being made, so that a genomic region can be joined in different combinations.")
parser$add_argument("-p", "--process_num", dest="process_num", action="store", type="integer", default=4, help="Number of processes to run with")
parser$add_argument("-f", "--force", dest="force", action="store_true", default=FALSE, help="Force overwrite")
user_args <- parser$parse_args()

# logging -----------------------------------------------------------------
library(log4r, quietly = T, warn.conflicts = F)
library(data.table, quietly = T)
library(dplyr, quietly = T, warn.conflicts = F)

my_layout <- function(level, ...) {
  paste0("[",format(Sys.time()), "]\t", user_args$group_name, "\t", level, "\t", Sys.getpid(), "\t", ..., "\n", collapse = "")
}

logger <- logger(threshold = "DEBUG", appenders = file_appender(file.path(dirname(user_args$input_no_snps), format(Sys.time(), "PoolEditingSites.%Y-%m-%dT%T.log")),
                                                                layout = my_layout))

debug(logger = logger, paste("User arguments:", user_args))

# External Scripts --------------------------------------------------------

# General -----------------------------------------------------------------

### Output directory
outdir = user_args$output_dir
force_run = user_args$force
num_processes = user_args$process_num

# non-SNP sites
nt_count_per_ES_pooled_noSNP = user_args$input_no_snps

# SNP sites
nt_count_per_ES_pooled_onlySNP = user_args$input_only_snps

# group size
current_group_size = user_args$current_group_size

debug(logger, paste("File without SNP:", nt_count_per_ES_pooled_noSNP))
debug(logger, paste("File with only SNP:", nt_count_per_ES_pooled_onlySNP))
debug(logger, paste("Group size:", current_group_size))

# create directory if does not exist
if(!dir.exists(outdir)) {
  dir.create(outdir, recursive = T)
}

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

# add a TRUE\FALSE column to count and filter by
add_filters <- function(df, col_name, min_vec) {
  for (m in min_vec) {
    df = mutate(.data = df, "NumNonSNPASitesInRegionWith_{{col_name}}AtLeast{{m}}" := if_else({{col_name}} >= m, 
                                                                                              true = TRUE,
                                                                                              false = FALSE, 
                                                                                              missing = FALSE))
  }
  return(df)
}

add_filters_exclusive <- function(df, col_name, min_vec) {
  for (m in min_vec) {
    df = mutate(.data = df, "NumNonSNPASitesInRegionWith_{{col_name}}Over{{m}}" := if_else({{col_name}} > m, 
                                                                                              true = TRUE,
                                                                                              false = FALSE, 
                                                                                              missing = FALSE))
  }
  return(df)
}

## add combinations of conditions - coverage count and minimal editing
add_combinations <- function(df, editing_col, cov_col, min_editing_vec, min_cov_vec) {
  for (min_editing in min_editing_vec) {
    for (min_G_cov in min_cov_vec) {
      df = mutate(.data = df, "NumNonSNPASitesInRegionWith_{{editing_col}}AtLeast{{min_editing}}_Num{{cov_col}}AtLeast{{min_G_cov}}" := 
                    if_else({{editing_col}} >= min_editing & {{cov_col}} >= min_G_cov, 
                            true = TRUE,
                            false = FALSE, 
                            missing = FALSE))
    }
  }
  return(df)
}

add_combinations_exclusive_editing <- function(df, editing_col, cov_col, min_editing_vec, min_cov_vec) {
  for (min_editing in min_editing_vec) {
    for (min_G_cov in min_cov_vec) {
      df = mutate(.data = df, "NumNonSNPASitesInRegionWith_{{editing_col}}Over{{min_editing}}_Num{{cov_col}}AtLeast{{min_G_cov}}" := 
                    if_else({{editing_col}} > min_editing & {{cov_col}} >= min_G_cov, 
                            true = TRUE,
                            false = FALSE, 
                            missing = FALSE))
    }
  }
  return(df)
}

## add combinations of conditions - coverage count and minimal editing
add_single_combination <- function(df, editing_col, cov_col, min_editing, max_editing, min_cov) {
  df = mutate(.data = df, "NumNonSNPASitesInRegionWith_{{editing_col}}AtLeast{{min_editing}}_{{editing_col}}Under{{max_editing}}_Num{{cov_col}}AtLeast{{min_cov}}" := 
                if_else({{editing_col}} >= min_editing & {{editing_col}} < max_editing & {{cov_col}} >= min_cov, 
                        true = TRUE,
                        false = FALSE,
                        missing = FALSE))
}

add_single_combination_exculsive_editing <- function(df, editing_col, cov_col, min_editing, max_editing, min_cov) {
  df = mutate(.data = df, "NumNonSNPASitesInRegionWith_{{editing_col}}Over{{min_editing}}_{{editing_col}}Under{{max_editing}}_Num{{cov_col}}AtLeast{{min_cov}}" := 
                if_else({{editing_col}} > min_editing & {{editing_col}} < max_editing & {{cov_col}} >= min_cov, 
                        true = TRUE,
                        false = FALSE,
                        missing = FALSE))
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

# PROCESS =========================================================================================
# *************************************************************************************************


# process -------------------------------------------------------------
info(logger, "Processing file")

process_and_complement_files <- function(f) {
  ## plus -----
  # process files as is - assume all regions are on + strand and A2G is the wanted MM
  
  # if exists and should not force overwrite - skip
  if (file.exists(gsub(x = f, pattern = ".csv", replacement = ".processed.A2G.csv", fixed = T)) & !force_run) {
    warn(logger, paste("File already processed on plus strand, delete or force: ", gsub(x = f, pattern = ".csv", replacement = ".processed.A2G.csv", fixed = T)))
  } else {
    debug(logger, paste("Processing", f))
    fread(f, stringsAsFactors = F, header = T) %>% 
      # use as-is, no need to process reference or count column names
      process %>%
      fwrite(gsub(x = f, pattern = ".csv", replacement = ".processed.A2G.csv", fixed = T),
             quote = F, row.names = F, scipen = 999)
    # inform to ensure writing is finished
    debug(logger, paste("Finished writing", gsub(x = f, pattern = ".csv", replacement = ".processed.A2G.csv", fixed = T)))
    gc()
  }
  
  # minus -----
  # swap files to match complementary strand
  # then process new files - assume all regions are on - strand and T2C is the wanted MM
  if (file.exists(gsub(x = f, pattern = ".csv", replacement = ".processed.T2C.csv", fixed = T)) & !force_run) {
    warn(logger, paste("File already processed on minus strand, delete or force: ", gsub(x = f, pattern = ".csv", replacement = ".processed.T2C.csv", fixed = T)))
  } else {
    debug(logger, paste("Processing and complementing", f))
    fread(f, stringsAsFactors = F, header = T) %>% 
      # complement - take T reference and swap bases
      complement %>%
      # process as if plus after complementation
      process %>%
      fwrite(gsub(x = f, pattern = ".csv", replacement = ".processed.T2C.csv", fixed = T),
             quote = F, row.names = F, scipen = 999)
    # inform to ensure writing is finished
    debug(logger, paste("Finished writing", gsub(x = f, pattern = ".csv", replacement = ".processed.T2C.csv", fixed = T)))
    gc()
  }
}


# combine and pool --------------------------------------------------------

# combine\process the files
read_files <- function(f) {
  f = f %>%
    # read
    fread(stringsAsFactors = F, header = T) %>%
    mutate(Group = user_args$group_name)
  
  # if needed, change region ID
  if(user_args$region_id_file != "") {
    debug(logger, paste("Loading new region ID file", user_args$region_id_file))
    new_region_id = fread(user_args$region_id_file, stringsAsFactors = F, col.names = c("Region", "NewRegionID"))
    debug(logger, paste("Joining new region ID file", user_args$region_id_file))
    # info(logger, paste("Joining", unique(f$Region), "exisiting regions with", 
    #                    unique(new_region_id)), "new regions")
    f = f %>%
      # fix case of empty file
      mutate(Region = as.character(Region)) %>%
      # add new ids
      inner_join(new_region_id, by = "Region") %>%
      # remove old ones
      select(-Region) %>%
      # rename
      rename(Region = NewRegionID)
    debug(logger, paste("Finished joining new region ID file", user_args$region_id_file))
    # info(logger, paste("Analyzing", unique(f$Region), "regions left"))
  }
  
  return(f)
}


# function to pool no SNP (main info)
pool_no_snp <- function(HE_nt_count_per_ES_pooled_noSNP, chr_region, current_group_size, strand) {
  debug(logger, paste(strand, "\tPooling all bases in", chr_region))
  
  HE_nt_count_per_ES_pooled_noSNP %>%
    filter(PositionChr == chr_region) %>%
    group_by(Region, Reference, Group) %>%
    # group only by region, treat all nts together (each will only affect the relevant indexes)
    # group_by(Region) %>%
    summarise(TotalIndexedSitesInRegion = n(),
              # sum counts and coverage
              across(c(A, C, G, `T`, TotalCoverage), sum)) %>%
    # add editing index
    process() %>%
    rename_with(.cols = ends_with("EditingIndex"), ~paste0(.x, "NoSNP"))  %>%
    # add number of samples per group
    # add coverage per ES and per ES + sample
    mutate(TotalSamplesPerGroup = current_group_size,
           AvgAGCoveragePerSite = AGCoveragePerSite / TotalSamplesPerGroup,
           AvgTotalCoveragePerSite = TotalCoverage / TotalSamplesPerGroup,
           MeanCoveragePerSite = TotalCoverage / TotalIndexedSitesInRegion,
           AvgMeanCoveragePerSiteNoSNP = MeanCoveragePerSite / TotalSamplesPerGroup) %>%
    relocate(TotalSamplesPerGroup, .after = Group) %>% 
    # create one row per region
    tidyr::pivot_wider(names_from = Reference,
                       values_from = TotalIndexedSitesInRegion:AvgMeanCoveragePerSiteNoSNP,
                       names_glue = "{.value}_ReferenceBase{Reference}") %>%
    # remove the all-NA columns (for example A2G editing index for C reference base)
    # remove meaningless rows (such as AG coverage for non-A reference)
    select(where(~any(!is.na(.x))), -matches("AGCoveragePerSite_ReferenceBase[CGT]")) %>%
    # fix editing index columns
    rename_with(.cols = contains("EditingIndex"), ~gsub(.x, pattern = "_ReferenceBase[ACGT]", replacement = "")) %>% 
    # add coverage per ES and per ES + sample
    mutate(MeanAGCoveragePerSiteNoSNP_ReferenceBaseA = AGCoveragePerSite_ReferenceBaseA / TotalIndexedSitesInRegion_ReferenceBaseA,
           AvgMeanAGCoveragePerSiteNoSNP_ReferenceBaseA = MeanAGCoveragePerSiteNoSNP_ReferenceBaseA / TotalSamplesPerGroup,
           MeanCoveragePerSiteNoSNP = sum(c_across(matches("TotalCoverage_ReferenceBase[ACGT]"))) / sum(c_across(matches("TotalIndexedSitesInRegion_ReferenceBase[ACGT]"))),
           AvgMeanCoveragePerSiteNoSNP = MeanCoveragePerSiteNoSNP / TotalSamplesPerGroup) %>%
    rename_with(.cols = contains("AGCoverage"), ~gsub(.x, pattern = "_ReferenceBaseA", replacement = "")) %>%
    # reorder
    select(Region, contains("AGCoverage"), contains("EditingIndex"), starts_with("TotalIndexedSitesInRegion"),
           contains("Coverage"), everything(), ends_with("_ReferenceBaseN"))
}


# function to pool only A of no SNP
pool_no_snp_A <- function(HE_nt_count_per_ES_pooled_noSNP, HE_nt_count_per_ES_pooled_noSNP_aggregated, chr_region, current_group_size, strand) {
  debug(logger, paste(strand, "\tCounting site pooling for A bases only without SNP in", chr_region))
  
  HE_nt_count_per_ES_pooled_noSNP %>%
    # filter(Reference == "A") %>%
    filter(Reference == "A", PositionChr == chr_region) %>%
    # add number of samples per group
    # add coverage per ES and per ES + sample
    mutate(TotalSamplesPerGroup = current_group_size,
           AvgAGCoveragePerSite = AGCoveragePerSite / TotalSamplesPerGroup) %>%
    add_combinations(editing_col = A2GEditingIndex,
                     cov_col = G,
                     min_editing_vec = c(1, seq(5, 30, 5)),
                     min_cov_vec = c(2)) %>%
    add_combinations_exclusive_editing(editing_col = A2GEditingIndex,
                                       cov_col = G,
                                       min_editing_vec = c(0),
                                       min_cov_vec = c(2)) %>%
    add_single_combination_exculsive_editing(editing_col = A2GEditingIndex, 
                                             cov_col = G, 
                                             min_editing = 0, 
                                             max_editing = 1, 
                                             min_cov = 2) %>%
    add_single_combination(editing_col = A2GEditingIndex, 
                           cov_col = G, 
                           min_editing = 5, 
                           max_editing = 10, 
                           min_cov = 2) %>%
    add_single_combination(editing_col = A2GEditingIndex, 
                           cov_col = G, 
                           min_editing = 1, 
                           max_editing = 5, 
                           min_cov = 2) %>%
    add_filters(col_name = A2GEditingIndex, 
                min_vec = c(1, 5, 10, 25, 50, 75, 90)) %>%
    add_filters(col_name = AvgAGCoveragePerSite, 
                min_vec = c(0.25, 0.5, 1, 5, 10, 25)) %>%
    add_filters_exclusive(col_name = A2GEditingIndex, 
                          min_vec = c(0)) %>%
    # add global index
    inner_join(HE_nt_count_per_ES_pooled_noSNP_aggregated %>%
                 select(Region, A2GEditingIndexNoSNP) %>%
                 rename(GlobalA2GEditingIndexNoSNP = A2GEditingIndexNoSNP)) %>%
    # add for each site whether it is greater than global index
    mutate("NumNonSNPASitesInRegionWith_A2GEditingIndexAtLeastGlobalA2GEditingIndexNoSNP_NumGAtLeast2" =
             if_else(A2GEditingIndex >= GlobalA2GEditingIndexNoSNP & G >= 2, 
                     true = TRUE,
                     false = FALSE, 
                     missing = FALSE)) %>%
    group_by(Region, Group) %>%
    summarise(TotalIndexedSitesInRegion_ReferenceBaseA = n(),
              across(starts_with("NumNonSNPASitesInRegionWith"), sum),
              # add percentage
              across(starts_with("NumNonSNPASitesInRegion"), ~.x*100/TotalIndexedSitesInRegion_ReferenceBaseA, .names = "%{.col}")) %>%
    # fix names
    rename_with(.cols = starts_with("%NumNonSNPASitesInRegion"),
                ~stringr::str_replace(., "Num", replacement = ""))
  
}


# function to pool only SNP
pool_only_snp <- function(HE_nt_count_per_ES_pooled_onlySNP, chr_region, strand) {
  debug(logger, paste(strand, "\tCounting site pooling for SNP bases only in", chr_region))
	
  HE_nt_count_per_ES_pooled_onlySNP %>%
    filter(PositionChr == chr_region) %>%
    group_by(Region, Reference, Group) %>%
    summarise(TotalSNPSitesInRegion_ReferenceBase = n()) %>%
    tidyr::pivot_wider(names_from = Reference,
                       values_from = TotalSNPSitesInRegion_ReferenceBase,
                       names_glue = "{.value}{Reference}", 
                       values_fill = 0) %>%
    mutate(TotalSNPSitesInRegion = sum(c_across(starts_with("TotalSNPSitesInRegion_ReferenceBase")), na.rm = T))
}


# one function to do it all
pool_by_chr <- function(nt_count_per_ES_pooled_noSNP, nt_count_per_ES_pooled_onlySNP, outdir, group_name, current_group_size, force_run = FALSE, suffix = "", strand) {
  # make sure this is not in vain
  if (file.exists(file.path(outdir, paste0("EditedRegionsSitesInformation.pooled", suffix, ".csv"))) & !force_run) {
    warn(logger, paste("Final file", file.path(outdir, paste0("EditedRegionsSitesInformation.pooled", suffix, ".csv")), "exists, delete or force. Skipping..."))
    return("Exit")
  }
  
  
  # READ =========================================================================================
  # *************************************************************************************************
  info(logger, paste0(strand, "\tLoading no-SNP file ", nt_count_per_ES_pooled_noSNP))# read files
  HE_nt_count_per_ES_pooled_noSNP = nt_count_per_ES_pooled_noSNP %>%
    read_files
  
  # POOL REGION (NO SNP) ============================================================================
  # *************************************************************************************************
  # pool all ES per region without SNPs

  # get chr -----------------------------------------------------------------
  # get all chromosomes
  chromosomes_vec = HE_nt_count_per_ES_pooled_noSNP %>%
    pull(PositionChr) %>%
    unique()
  info(logger, paste(strand, "\tFound the following chromosomes in no-SNP file:", paste(chromosomes_vec, collapse = "\n")))
  
  
  # pool information for all bases  -----------------------------------------
  if (!file.exists(file.path(outdir, paste0("EditedRegionsSites.pooled", suffix, ".csv"))) | force_run) {
    info(logger, paste0(strand, "\tStarting pooling for all bases"))
    
    # call
    HE_nt_count_per_ES_pooled_noSNP_aggregated = furrr::future_map_dfr(.x = chromosomes_vec,
                                                                       .f = ~ pool_no_snp(HE_nt_count_per_ES_pooled_noSNP = HE_nt_count_per_ES_pooled_noSNP,
                                                                                          current_group_size = current_group_size,
                                                                                          chr_region = .x, strand = strand))
    # HE_nt_count_per_ES_pooled_noSNP_aggregated = HE_nt_count_per_ES_pooled_noSNP %>%
    #   pool_no_snp
    fwrite(HE_nt_count_per_ES_pooled_noSNP_aggregated, file.path(outdir, paste0("EditedRegionsSites.pooled", suffix, ".csv")), quote = F, row.names = F, scipen = 999)
    # inform to ensure writing is finished
    debug(logger, paste(strand, "\tFinished writing", file.path(outdir, paste0("EditedRegionsSites.pooled", suffix, ".csv"))))
  } else {
    #load
    info(logger, paste0(strand, "\tLoading pooling for all bases"))
    HE_nt_count_per_ES_pooled_noSNP_aggregated = fread(file.path(outdir, paste0("EditedRegionsSites.pooled", suffix, ".csv")),
                                                       stringsAsFactors = F, header = T)
  }
  
  # pooling count of sites for A nts only -----------------------------------
  info(logger, paste0(strand, "\tStarting site count pooling for A bases only without SNP"))
  
  
  HE_nt_count_per_ES_pooled_noSNP_aggregated_A = furrr::future_map_dfr(.x = chromosomes_vec,
                                                                     .f = ~ pool_no_snp_A(HE_nt_count_per_ES_pooled_noSNP_aggregated = HE_nt_count_per_ES_pooled_noSNP_aggregated,# %>% filter(PositionChr == .x),
                                                                                          HE_nt_count_per_ES_pooled_noSNP = HE_nt_count_per_ES_pooled_noSNP,
                                                                                        current_group_size = current_group_size,
                                                                                        chr_region = .x, strand = strand))

  debug(logger, paste0(strand, "\tCleanup of information without SNP"))
  rm(HE_nt_count_per_ES_pooled_noSNP)
  gc()

  # pooling count of SNP -------------------------------------------------------
  info(logger, paste0(strand, "\tLoading SNP only file"))
  HE_nt_count_per_ES_pooled_onlySNP = nt_count_per_ES_pooled_onlySNP %>%
    read_files
  
  info(logger, paste0(strand, "\tStarting site count pooling for SNP bases only"))
  HE_nt_count_per_ES_pooled_onlySNP_aggregated = furrr::future_map_dfr(.x = chromosomes_vec,
                                                                       .f = ~ pool_only_snp(HE_nt_count_per_ES_pooled_onlySNP = HE_nt_count_per_ES_pooled_onlySNP,#  %>% filter(PositionChr == .x),
                                                                                            chr_region = .x, strand = strand))
  # HE_nt_count_per_ES_pooled_onlySNP_aggregated = HE_nt_count_per_ES_pooled_onlySNP %>%
  #   pool_only_snp
  
  debug(logger, paste0(strand, "\tCleanup of information of only SNP"))
  rm(HE_nt_count_per_ES_pooled_onlySNP)
  gc()
  
  
  # join and write output ---------------------------------------------------
  info(logger, paste0(strand, "\tJoining tables"))
  
  inner_join(HE_nt_count_per_ES_pooled_noSNP_aggregated,
             HE_nt_count_per_ES_pooled_noSNP_aggregated_A) %>%
    left_join(HE_nt_count_per_ES_pooled_onlySNP_aggregated) %>%
    tidyr::replace_na(list(TotalSNPSitesInRegion_ReferenceBaseA = 0, 
                           TotalSNPSitesInRegion = 0)) %>%
    fwrite(file.path(outdir, paste0("EditedRegionsSitesInformation.pooled", suffix, ".csv")), quote = F, row.names = F, scipen = 999)
  # inform to ensure writing is finished
  debug(logger, paste(strand, "\tFinished writing", file.path(outdir, paste0("EditedRegionsSitesInformation.pooled", suffix, ".csv"))))
  
}

pool_strand <- function(strand) {
  info(logger, paste(strand, "\tStarting pooling site information"))
  # fix file paths to match suffix
  pool_by_chr(nt_count_per_ES_pooled_noSNP = nt_count_per_ES_pooled_noSNP %>%
                     gsub(pattern = ".csv", replacement = paste0(".processed.", strand,".csv"), fixed = T),
                   nt_count_per_ES_pooled_onlySNP = nt_count_per_ES_pooled_onlySNP %>%
                     gsub(pattern = ".csv", replacement = paste0(".processed.", strand, ".csv"), fixed = T), 
                   # create suffix: group.A2G
                   suffix = paste0(".", user_args$group_name, ".", strand), 
                   group_name = user_args$group_name, current_group_size = current_group_size,
                   outdir = outdir, force_run = force_run, strand = strand)
  info(logger, paste0(strand, "\tFinished pooling site information"))
}

# main  -------------------------------------------------------------------
# set processing plan to multicore and set number of cores
options(future.globals.maxSize= 70000 * 1024 ^ 2)
future::plan(future::multicore, workers = num_processes)


# call processing
furrr::future_walk(.x = c(nt_count_per_ES_pooled_noSNP, nt_count_per_ES_pooled_onlySNP),
                   .f = ~ process_and_complement_files(f = .x))

# plus
pool_strand(strand = "A2G")
# minus
pool_strand(strand = "T2C")


info(logger, paste("Done pooling site information for", user_args$group_name))


# **************************************************************************

