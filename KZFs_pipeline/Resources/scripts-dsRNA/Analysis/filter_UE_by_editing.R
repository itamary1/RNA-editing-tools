
# Setup -------------------------------------------------------------------
rm(list = ls(all = TRUE))
library(data.table)
library(dplyr)

HE_EI_perRegionPerSample_pooled=fread("/private9/Projects/dsRNAProject/Summary/RNAEditingIndexOnHERegion/A2G/HyperEditedRegionsPerRegionPerSample.AllInfo.short.csv", 
                                      stringsAsFactors = F)


# filter by basic conditions
HE_EI_perRegionPerSample_pooled_noAlt = HE_EI_perRegionPerSample_pooled %>% 
  filter(!grepl(HERegion, pattern = "_")) %>%  # no alt
  select(HERegion, Length, Strand, Overlap, A2GEditingIndex, A2GEditingIndexNoSNP, 
         NonA2GEditingIndexSumNoSNPOverlapConsidered, EditingPcntNoSNPClass,
         starts_with("NumNonSNPASites") & ends_with("_NumGAtLeast2"), 
         AvgMeanCoveragePerRegion, starts_with("Total"),
         contains("Repeat"), contains("Genomic"), contains("Count"),
         ends_with("EditingIndex"), everything()) 


rm(HE_EI_perRegionPerSample_pooled)

# filter 1 ------------------------------------------------------------------

# take regions with at least 3 sites with 10 % editing and 
suffix = "_min5A2GEditingIndexNoSNP_max0.5NonA2GEditingIndexSumNoSNPOverlapConsidered_min3NumASitesInRegionWithA2GEditingIndexAtLeast5AndNumGAtLeast2"
HE_EI_perRegionPerSample_pooled_noAlt_min5A2GEditingIndexNoSNP_max0.5NonA2GEditingIndexSumNoSNPOverlapConsidered_min3NumASitesInRegionWithA2GEditingIndexAtLeast5AndNumGAtLeast2 = 
  HE_EI_perRegionPerSample_pooled_noAlt %>% 
  filter(A2GEditingIndexNoSNP >= 5, 
         NonA2GEditingIndexSumNoSNPOverlapConsidered <= 0.5, 
         NumNonSNPASitesInRegionWith_A2GEditingIndexAtLeast5_NumGAtLeast2 >= 3) %>% 
  arrange(-A2GEditingIndexNoSNP) 

fwrite(HE_EI_perRegionPerSample_pooled_noAlt_min5A2GEditingIndexNoSNP_max0.5NonA2GEditingIndexSumNoSNPOverlapConsidered_min3NumASitesInRegionWithA2GEditingIndexAtLeast5AndNumGAtLeast2, 
       file = paste0("/private9/Projects/dsRNAProject/Summary/RNAEditingIndexOnHERegion/A2G/HyperEditedRegionsPerRegionPerSample", suffix, ".AllInfo.short.csv"), 
       quote = F, row.names = F, scipen = 999)
fwrite(HE_EI_perRegionPerSample_pooled_noAlt_min5A2GEditingIndexNoSNP_max0.5NonA2GEditingIndexSumNoSNPOverlapConsidered_min3NumASitesInRegionWithA2GEditingIndexAtLeast5AndNumGAtLeast2 %>% 
         select(HERegion),
       file = paste0("/private9/Projects/dsRNAProject/Summary/RNAEditingIndexOnHERegion/A2G/Regions", suffix, ".pooled.txt"), 
       quote = F, row.names = F, scipen = 999,
       sep = "\t", col.names = F)
fwrite(HE_EI_perRegionPerSample_pooled_noAlt_min5A2GEditingIndexNoSNP_max0.5NonA2GEditingIndexSumNoSNPOverlapConsidered_min3NumASitesInRegionWithA2GEditingIndexAtLeast5AndNumGAtLeast2 %>% 
         select(HERegion, Length, Overlap, Strand) %>% 
         tidyr::separate(HERegion, into = c("Chr", "Start", "End"), sep = "[:-]") %>%
         mutate(Strand = recode(Strand, Plus = "+", Minus = "-")) %>%
         select(Chr, Start, End, Length, Overlap, Strand),
       file = paste0("/private9/Projects/dsRNAProject/Summary/RNAEditingIndexOnHERegion/A2G/Regions", suffix, ".pooled.bed"), 
       quote = F, row.names = F, scipen = 999,
       sep = "\t", col.names = F)



# > analyze -----------------------------------------------------------------

HE_EI_perRegionPerSample_pooled_noAlt_min5A2GEditingIndexNoSNP_max0.5NonA2GEditingIndexSumNoSNPOverlapConsidered_min3NumASitesInRegionWithA2GEditingIndexAtLeast5AndNumGAtLeast2 %>% 
  select(HERegion:`TissueCount_Edited%`, NonA2GEditingIndexSumNoSNPOverlapConsidered, NumNonSNPASitesInRegionWith_A2GEditingIndexAtLeast10_NumGAtLeast2, 
         ends_with("EditingIndex"), ends_with("EditingIndexNoSNP")) %>%
  View("Pass filters")
HE_EI_perRegionPerSample_pooled_noAlt %>% 
  select(HERegion:`TissueCount_Edited%`, NonA2GEditingIndexSumNoSNPOverlapConsidered, NumNonSNPASitesInRegionWith_A2GEditingIndexAtLeast10_NumGAtLeast2, 
         ends_with("EditingIndex"), ends_with("EditingIndexNoSNP")) %>%
  View("All")
HE_EI_perRegionPerSample_pooled_noAlt %>% 
  filter(!HERegion %in% HE_EI_perRegionPerSample_pooled_noAlt_min5A2GEditingIndexNoSNP_max0.5NonA2GEditingIndexSumNoSNPOverlapConsidered_min3NumASitesInRegionWithA2GEditingIndexAtLeast5AndNumGAtLeast2$HERegion) %>% 
  select(HERegion:`TissueCount_Edited%`, NonA2GEditingIndexSumNoSNPOverlapConsidered, NumNonSNPASitesInRegionWith_A2GEditingIndexAtLeast10_NumGAtLeast2, 
         ends_with("EditingIndex"), ends_with("EditingIndexNoSNP")) %>%
  View("Failed to pass filters")




# Top 30 good
HE_EI_perRegionPerSample_pooled_noAlt_min5A2GEditingIndexNoSNP_max0.5NonA2GEditingIndexSumNoSNPOverlapConsidered_min3NumASitesInRegionWithA2GEditingIndexAtLeast5AndNumGAtLeast2 %>%  
  select(HERegion:NonA2GEditingIndexSumNoSNPOverlapConsidered, NumNonSNPASitesInRegionWith_A2GEditingIndexAtLeast10_NumGAtLeast2, TotalSNPSitesInRegion) %>%
  arrange(-A2GEditingIndexNoSNP) %>%
  head(n=30) %>% 
  View("Top 30 Pass filters")

# Top 30 bad
HE_EI_perRegionPerSample_pooled_noAlt %>% 
  filter(!HERegion %in% HE_EI_perRegionPerSample_pooled_noAlt_min5A2GEditingIndexNoSNP_max0.5NonA2GEditingIndexSumNoSNPOverlapConsidered_min3NumASitesInRegionWithA2GEditingIndexAtLeast5AndNumGAtLeast2$HERegion) %>% 
  select(HERegion:NonA2GEditingIndexSumNoSNPOverlapConsidered, NumNonSNPASitesInRegionWith_A2GEditingIndexAtLeast10_NumGAtLeast2, TotalSNPSitesInRegion) %>%
  arrange(-A2GEditingIndexNoSNP) %>% 
  head(n=30) %>% 
  View("Top 30 Fail filters")


# *********************************************************************************************************************************

