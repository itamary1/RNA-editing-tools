###########################################################################################################################
All this dir and subdir content was written by Roni Cohen for her own usage. I copied it as is, and add RID_creator + change  the pathes
###########################################################################################################################

# scripts-dsRNA
Pipelines and scripts for preprocess, pooling and analysis of region data after run of [RNA editing index](https://github.com/a2iEditing/RNAEditingIndexer).

## 	Calculate RNA editing index
Self-done, including `--per_region`, `--per_sample`, `--keep_cmpileup` parameters.

## Processing
### Pool information per region by groups
*For each group, pool all region information*     

* Script to aggregate\split regions information by chromosomes: scripts-dsRNA/Processing/Editing/split_StrandDerivingCountsPerRegion_file_by_chromosomes.R    
* Script to pool information by group: scripts-dsRNA/Processing/Editing/Region_perRegionPerSample_GTExSubset_preprocess.withSignal.R     
* **Pipeline for basic information post-processing scripts: scripts-dsRNA/Pipelines/post_editing_index_calculation_analysis.sh**         

*Algorithm*
* Split each sample into one united file per chromosome
    * Editing index information
        * First change the name of the top-level file (as it will be caught by the suffix) then run uniting script
* Pooling script: 
    * Editing index information
        * pool each chromosome per group
        * Combine processed files per group
    * Indexed sites information
		* Pool each group
(sum is commutative, without sample information due to memory limits)
		
### Extract number of canonical and mismatch events per site in each region
*For all bases, remove SNP positions*      

* Script to aggregate nt count per ES (sum cmpileup for all positions): scripts-dsRNA/Processing/Cmpileup/sum_cmpileup_allMM.py       
* Script to split into SNP and non-SNP sites:	scripts-dsRNA/Processing/Cmpileup/summed_cmpileup_filterSNP.sh     
* Script to aggregate ES information (cmpileup script): scripts-dsRNA/Processing/Cmpileup/all_ES_in_regions_analysis.R    
* **Pipeline for cmpileup post-processing scripts: scripts-dsRNA/Pipelines/post_RNA_editing_index_cmpileup_analysis.sh**     

## Post-processing:
*Join all information and filter regions*      

* **Script to aggregate region information (standalone): scripts-dsRNA/Analysis/all_UE_regions_analysis.R**       
* **Script to filter regions: scripts-dsRNA/Analysis/filter_UE_by_editing.R**        

*Algorithm*
* Combine all information for each region:
    * Editing index info
    * Count of edited samples and group for each region + group information
    * Percentage of total samples/group calculation
    * Category classification for editing and coverage
    * Repeat and genomic location annotation
    * Indexed sites pooled information (below)
* Indexed sites info
    * Editing & coverage
        * Non-SNP
        * SNP
    * Number of A sites passing thresholds within region
        * Non-SNP
        * SNP