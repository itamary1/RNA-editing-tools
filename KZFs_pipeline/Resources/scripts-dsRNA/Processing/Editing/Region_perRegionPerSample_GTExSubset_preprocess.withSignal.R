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
parser$add_argument("-c", "--columns_all", dest="group_file", action="store", required=TRUE, help="Column names for first file type")
parser$add_argument("-c1", "--columns_other_MM1", dest="group_file", action="store", required=TRUE, help="column names for second file type")
parser$add_argument("-c2", "--columns_other_MM2", dest="group_file", action="store", required=TRUE, help="CSV with sample to group. Will pool samples according to given groups, sample count per group will be calculated according to this file. Format: Sample,Group")
parser$add_argument("-g", "--group_file", dest="group_file", action="store", required=TRUE, help="CSV with sample to group. Will pool samples according to given groups, sample count per group will be calculated according to this file. Format: Sample,Group")
user_args <- parser$parse_args()

print(paste("User arguments:", user_args))

# Sourcing ----------------------------------------------------------------

# General -----------------------------------------------------------------

### Per region per sample, Alu index (via Hillel's pipeline file, summary files: "EditingIndex.csv")
EI_perRegionPerSample = list.files(user_args$input_directory, 
                                      pattern = "StrandDerivingCountsPerRegion.chr.*.preprocessed.strand.csv",
                                      full.names = T, recursive = F)
EI_perRegionPerSample_1 = list.files(user_args$input_directory, 
                                        pattern = "StrandDerivingCountsPerRegion.chr.*.preprocessed.otherMM1.strand.csv",
                                        full.names = T, recursive = F)
EI_perRegionPerSample_2 = list.files(user_args$input_directory, 
                                        pattern = "StrandDerivingCountsPerRegion.chr.*.preprocessed.otherMM2.strand.csv",
                                        full.names = T, recursive = F)


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
library(data.table, quietly = T)
library(dtplyr, quietly = T)
library(dplyr, warn.conflicts = FALSE, quietly = T)
sra = read.csv(file=user_args$group_file, header=T, check.names = F, stringsAsFactors=FALSE,
               col.names = c("Sample", "Group")) %>%
  # add group count
  group_by(Group) %>% 
  add_tally(name = "TotalSamplesPerGroup")  

# colors -------------------------------------------------------------------
# GTEx_colors <- read.delim(GTEx_colors, header = T, check.names = F, stringsAsFactors = F) 



# column names ------------------------------------------------------------
column_names = c("Chromosome", "Sample", "AluElement", "A2G", "Length", 
                 "SenseGeneCommonName", "SenseGeneRefSeqID", "AntisenseGenomicPosition", "SenseGenomicPosition", "AntisenseGeneRefSeqID", "AntisenseGeneCommonName", 
                 "A2GEditingIndex", "C2TEditingIndex", "TotalCoverageAtAllPositions", "MeanCoveragePerRegion", 
                 "IndexedMismatchesOfA2G", "IndexedCanonicalOfA2G", "IndexedMismatchesOfC2T", "IndexedCanonicalOfC2T",
                 "NumOfIndexedMismatchesSitesOfA2G", "NumOfIndexedOverallSitesOfA2G", "NumOfIndexedMismatchesSitesOfC2T", "NumOfIndexedOverallSitesOfC2T")

column_names1 = c("Chromosome", "Sample", "AluElement", "A2G", "Length", "TotalCoverageAtAllPositions", 
                  "IndexedCanonicalOfC2A", "IndexedCanonicalOfG2C", "IndexedCanonicalOfA2G", 
                  "IndexedMismatchesOfC2G", "IndexedMismatchesOfC2A", "IndexedCanonicalOfG2T", 
                  "IndexedMismatchesOfG2C", "IndexedMismatchesOfG2T", "NumOfIndexedOverallSitesOfC2G", 
                  "NumOfIndexedOverallSitesOfC2A", "IndexedCanonicalOfC2G", "NumOfIndexedOverallSitesOfG2C",
                  "NumOfIndexedOverallSitesOfG2T", "NumOfIndexedOverallSitesOfA2G", "NumOfIndexedMismatchesSitesOfG2C", 
                  "NumOfIndexedMismatchesSitesOfG2T", "IndexedMismatchesOfA2G", "NumOfIndexedMismatchesSitesOfC2A", 
                  "NumOfIndexedMismatchesSitesOfC2G", "NumOfIndexedMismatchesSitesOfA2G", "IndexedCanonicalOfG2A", 
                  "NumOfIndexedOverallSitesOfG2A", "IndexedMismatchesOfG2A", "NumOfIndexedMismatchesSitesOfG2A")

column_names2 = c("Chromosome", "Sample", "AluElement", "A2G", "Length", "TotalCoverageAtAllPositions", 
                  "IndexedCanonicalOfA2G", "IndexedCanonicalOfA2C", "IndexedCanonicalOfA2T", "NumOfIndexedMismatchesSitesOfT2A", 
                  "NumOfIndexedMismatchesSitesOfT2G", "IndexedMismatchesOfA2C", "IndexedMismatchesOfA2T", 
                  "IndexedCanonicalOfT2C", "IndexedCanonicalOfT2A", "NumOfIndexedMismatchesSitesOfT2C", 
                  "NumOfIndexedOverallSitesOfA2G", "NumOfIndexedOverallSitesOfA2C", "NumOfIndexedOverallSitesOfA2T", 
                  "NumOfIndexedOverallSitesOfT2G", "NumOfIndexedOverallSitesOfT2C", "NumOfIndexedOverallSitesOfT2A", 
                  "IndexedMismatchesOfA2G", "IndexedMismatchesOfT2G", "IndexedMismatchesOfT2A", "IndexedMismatchesOfT2C", 
                  "NumOfIndexedMismatchesSitesOfA2T", "NumOfIndexedMismatchesSitesOfA2C", "NumOfIndexedMismatchesSitesOfA2G", "IndexedCanonicalOfT2G")


# grouping ----------------------------------------------------------------
group_by_vars=c("AluElement", "A2G", "Group")

# ANALYSIS =======================================================================================
# *************************************************************************************************

# pooling of core date ---------------------------------------------------------------
for (f in c(EI_perRegionPerSample)) {
  # get chromosome number
  chromosome = gsub(basename(f), pattern = "StrandDerivingCountsPerRegion.|.preprocessed.csv", replacement = "")
  print(paste(chromosome, f))
  
  EI_perRegionPerSample <- fread(file = f, header = F, stringsAsFactors = F,
                                             col.names =  column_names) %>%
    inner_join(sra) %>%
    as.data.table()
  # fwrite(EI_perRegionPerSample, file.path(out_stats, paste0("InvertedElementsPerRegionPerSample.withSignal.csv")), quote = F, row.names = F, scipen = 999)
  
  print(typeof(EI_perRegionPerSample))
  print(paste("Finished loading table", f))
  print(head(EI_perRegionPerSample))
  
  # change suffix
  output_suffix = paste0(chromosome, ".")
  print(paste("Beginning pooling analysis for", chromosome))
  
  ## pool -------------------------------------------
  # pool by element
  EI_perRegionPerSample_pooled = EI_perRegionPerSample %>%
    # lazy_dt() %>%
    group_by(across(all_of(group_by_vars))) %>%
    summarise(across(A2GEditingIndex:IndexedCanonicalOfC2T, sum),
              across(c(Length:AntisenseGeneCommonName, TotalSamplesPerGroup), first),
              across(matches("NumOfIndexed(Mismatches|Overall)SitesOf[ACGT]2[ACGT]"), list(Min=min, Max=max), .names = "{.fn}{.col}"),
              SampleCount = n()) %>%
    mutate(MeanCoveragePerRegion = TotalCoverageAtAllPositions / Length,
           A2GEditingIndex = 100 * IndexedMismatchesOfA2G /  (IndexedMismatchesOfA2G + IndexedCanonicalOfA2G),
           C2TEditingIndex = 100 * IndexedMismatchesOfC2T /  (IndexedMismatchesOfC2T + IndexedCanonicalOfC2T),
           AvgMeanCoveragePerRegionInGroup = MeanCoveragePerRegion / TotalSamplesPerGroup) %>%
    select(AluElement, Length, contains("EditingIndex"), AvgMeanCoveragePerRegionInGroup, MeanCoveragePerRegion, 
           SampleCount,
           SenseGeneCommonName, AntisenseGeneCommonName, everything()) %>%
    as.data.table()
  print(paste("Writing", f, "into", file.path(out_stats, paste0("EditingPerRegionPerSample.", output_suffix,"pooled.withSignal.csv"))))
  fwrite(EI_perRegionPerSample_pooled, file.path(out_stats, paste0("EditingPerRegionPerSample.", output_suffix,"pooled.withSignal.csv")), quote = F, row.names = F, scipen = 999)
  
  ## sample and Group information -------------------------------------------
  # pool by element
  EI_perRegionPerSample_pooled_count = EI_perRegionPerSample %>%
    # lazy_dt() %>%
    group_by(across(all_of(group_by_vars))) %>%
    summarise(SampleCount = n(),
              SampleCountEdited = sum(A2GEditingIndex > 0, na.rm = T)) %>%
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
  print(paste("Writing", f, "into", file.path(out_stats, paste0("CountDetailsEditingPerRegionPerSample.", output_suffix,"pooled.withSignal.csv"))))
  fwrite(EI_perRegionPerSample_pooled_count, file.path(out_stats, paste0("CountDetailsEditingPerRegionPerSample.", output_suffix,"pooled.withSignal.csv")), quote = F, row.names = F, scipen = 999)
}




# pooling otherMM1 --------------------------------------------------------
for (f in c(EI_perRegionPerSample_1)) {
  # get chromosome number
  chromosome = gsub(basename(f), pattern = "StrandDerivingCountsPerRegion.|.preprocessed.otherMM1.csv", replacement = "")
  # get part
  part = gsub(basename(f), pattern = "StrandDerivingCountsPerRegion.chr.*.preprocessed.|.csv", replacement = "")
  print(paste(chromosome, part))
  
  EI_perRegionPerSample <- fread(file = f, header = F, stringsAsFactors = F,
                                    col.names =  column_names1) %>%
    inner_join(sra) %>%
    as.data.table()
  # fwrite(EI_perRegionPerSample, file.path(out_stats, paste0("InvertedElementsPerRegionPerSample.",  part, ".csv")), quote = F, row.names = F, scipen = 999)
  
  print(typeof(EI_perRegionPerSample))
  print(paste("Finished loading table", f))
  print(head(EI_perRegionPerSample))
  
  output_suffix = paste0(chromosome, ".")
  print(paste("Beginning pooling analysis for", output_suffix, part))
  
  # pool by element
  EI_perRegionPerSample_pooled_1 = EI_perRegionPerSample %>%
    # lazy_dt() %>%
    group_by(across(all_of(group_by_vars))) %>%
    summarise(across(matches("IndexedCanonicalOf[ATCG]2[ATCG]|IndexedMismatchesOf[ATCG]2[ATCG]|TotalCoverageAtAllPositions"), sum),
              across(c(Length, TotalSamplesPerGroup), first),
              across(contains("Sites"), list(Min=min, Max=max), .names = "{.fn}{.col}"),
              SampleCount = n()) %>%
    mutate(MeanCoveragePerRegion = TotalCoverageAtAllPositions / Length,
           A2GEditingIndex = 100 * IndexedMismatchesOfA2G /  (IndexedMismatchesOfA2G + IndexedCanonicalOfA2G),
           C2GEditingIndex = 100 * IndexedMismatchesOfC2G /  (IndexedMismatchesOfC2G + IndexedCanonicalOfC2G),
           C2AEditingIndex = 100 * IndexedMismatchesOfC2A /  (IndexedMismatchesOfC2A + IndexedCanonicalOfC2A),
           G2CEditingIndex = 100 * IndexedMismatchesOfG2C /  (IndexedMismatchesOfG2C + IndexedCanonicalOfG2C),
           G2TEditingIndex = 100 * IndexedMismatchesOfG2T /  (IndexedMismatchesOfG2T + IndexedCanonicalOfG2T),
           G2AEditingIndex = 100 * IndexedMismatchesOfG2A /  (IndexedMismatchesOfG2A + IndexedCanonicalOfG2A),
           AvgMeanCoveragePerRegionInGroup = MeanCoveragePerRegion / TotalSamplesPerGroup) %>%
    select(AluElement, Length, contains("EditingIndex"), AvgMeanCoveragePerRegionInGroup, MeanCoveragePerRegion, 
           SampleCount, everything()) %>%
    as.data.table()
  print(paste("Writing", f, "into", file.path(out_stats, paste0("EditingPerRegionPerSample.", output_suffix,"pooled.",  part, ".csv"))))
  fwrite(EI_perRegionPerSample_pooled_1, file.path(out_stats, paste0("EditingPerRegionPerSample.", output_suffix,"pooled.",  part, ".csv")), quote = F, row.names = F, scipen = 999)
}


# pooling other MM2 -------------------------------------------------------
for (f in c(EI_perRegionPerSample_2)) {
  # get chromosome number
  chromosome = gsub(basename(f), pattern = "StrandDerivingCountsPerRegion.|.preprocessed.otherMM2.csv", replacement = "")
  # get part
  part = gsub(basename(f), pattern = "StrandDerivingCountsPerRegion.chr.*.preprocessed.|.csv", replacement = "")
  print(paste(chromosome, part))
  
  EI_perRegionPerSample <- fread(file = f, header = F, stringsAsFactors = F,
                                    col.names =  column_names2) %>%
    inner_join(sra) %>%
    as.data.table()
  # fwrite(EI_perRegionPerSample, file.path(out_stats, paste0("InvertedElementsPerRegionPerSample.",  part, ".csv")), quote = F, row.names = F, scipen = 999)
  
  print(typeof(EI_perRegionPerSample))
  print(paste("Finished loading table", f))
  print(head(EI_perRegionPerSample))
  
  output_suffix = paste0(chromosome, ".")
  print(paste("Beginning pooling analysis for", output_suffix, part))
  
  # pool by element and group
  EI_perRegionPerSample_pooled_2 = EI_perRegionPerSample %>%
    # lazy_dt() %>%
    group_by(across(all_of(group_by_vars))) %>%
    summarise(across(matches("IndexedCanonicalOf[ATCG]2[ATCG]|IndexedMismatchesOf[ATCG]2[ATCG]|TotalCoverageAtAllPositions"), sum),
              across(c(Length, TotalSamplesPerGroup), first),
              across(contains("Sites"), list(Min=min, Max=max), .names = "{.fn}{.col}"),
              SampleCount = n()) %>%
    mutate(MeanCoveragePerRegion = TotalCoverageAtAllPositions / Length,
           A2GEditingIndex = 100 * IndexedMismatchesOfA2G /  (IndexedMismatchesOfA2G + IndexedCanonicalOfA2G),
           A2CEditingIndex = 100 * IndexedMismatchesOfA2C /  (IndexedMismatchesOfA2C + IndexedCanonicalOfA2C),
           A2TEditingIndex = 100 * IndexedMismatchesOfA2T /  (IndexedMismatchesOfA2T + IndexedCanonicalOfA2T),
           T2GEditingIndex = 100 * IndexedMismatchesOfT2G /  (IndexedMismatchesOfT2G + IndexedCanonicalOfT2G),
           T2AEditingIndex = 100 * IndexedMismatchesOfT2A /  (IndexedMismatchesOfT2A + IndexedCanonicalOfT2A),
           T2CEditingIndex = 100 * IndexedMismatchesOfT2C /  (IndexedMismatchesOfT2C + IndexedCanonicalOfT2C),
           AvgMeanCoveragePerRegionInGroup = MeanCoveragePerRegion / TotalSamplesPerGroup) %>%
    select(AluElement, Length, contains("EditingIndex"), AvgMeanCoveragePerRegionInGroup, MeanCoveragePerRegion, 
           SampleCount, everything()) %>%
    as.data.table()
  print(paste("Writing", f, "into", file.path(out_stats, paste0("EditingPerRegionPerSample.", output_suffix,"pooled.",  part, ".csv"))))
  fwrite(EI_perRegionPerSample_pooled_2, file.path(out_stats, paste0("EditingPerRegionPerSample.", output_suffix,"pooled.",  part, ".csv")), quote = F, row.names = F, scipen = 999)
}


print("Done")

# # test - if this can be traded by an inner join of all (sort of for on "zip" like in python, or by path processing)
# inner_join(EI_perRegionPerSample, EI_perRegionPerSample_1) %>%
#   inner_join(EI_perRegionPerSample_2) %>%
#   fwrite(file.path(out_stats, paste0("EditingPerRegionPerSample.", output_suffix,"pooled.csv")), quote = F, row.names = F, scipen = 999)




# **********************************************************************************************************


