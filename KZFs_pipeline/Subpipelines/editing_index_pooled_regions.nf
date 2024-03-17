
params.help=false
params.have_base_config=false
if((!params.have_base_config) && (!params.help)){
 println "please attach a base config file"
 System.exit(1) 
}

if((!params.profile_selected) && (!params.help) && (!params.pooledEI_snps_file)){
println "please select a profile with -profile <proile_name> or set pooledEI_snps_file"
 System.exit(1) 
}

File file = new File(params.pooledEI_snps_file)
if((!file.exists()) && (!params.help)){
 println "file for pooledEI_snps_file not exist at $params.pooledEI_snps_file, please select a profile with -profile <proile_name> or set pooledEI_snps_file"
    System.exit(1)
}



// if((params.docker_run) && (!params.bams_dir) && (!params.help)){
//     println "you must set bams dir when runing with docker. exiting.."
//     System.exit(1)
// }

def helpMessage() {
    log.info '''\
            GET RNA EDITING INDEX ON POOLED REGIONS
            ==================================
            params to set:
              required:

                //the parent directory of all your bams
                bams_dir = ''

                //regions for cmpileup to run on - it can be bed6 which you can also use for the pool_bed6_file(see below)
                cmp_regions_bed = ''

                // RID_file or bed6 file is required (you dont need both)
                //region -> id file that map evey regoin to a group id (have no relation to the group_file  above) containig: 
                //1.two columns header(names doesnt matter)
                //2. in every line: a. Region (chr:start-end), b. group id (like gene-name)
                RID_file=''

                // you can also use my inner RID creator and supply bed6 file
                // column 4 must be the group id (gene name etc..) !!!
                pool_bed6_file=''

                recommended:
                pooledEI_result_dir = "$launchDir/PooledEditingIndex"
                bams_suffix='.Aligned.sortedByCoord.out.bam'


                other params:
                // group file - if you want the results to be for groups of samples (control/condition etc..)
                // if you want your result for every sample seperatly, you should leave this as empty string
                group_file_csv=''
                pool_CMP_max_parallell = 10
                which_python_39 = '/private/common/Software/anaconda/anaconda3/envs/python3/bin/python'
                which_pooling_script = "${params.special_pipelines_dir}/Resources/scripts-dsRNA/Pipelines/pool_RNA_editing_index.py"
                which_pooling_script_command = "$params.which_python_39 $params.which_pooling_script sites"
                which_createRID_command = "Rscript /home/alu/twerski/Scripts/Nextflow/Special_pipelines/Resources/scripts-dsRNA/RID_creator.R"
                profiles{
                hg38 {
                    params.profile_selected = "hg38"       
                    params.pooledEI_snps_file = "/private/common/Software/AEI/RNAEditingIndex1.1/RNAEditingIndexer/Resources/SNPs/HomoSapiens/ucscHg38CommonGenomicSNPs150.bed.gz"

                }
                mm10 {
                    params.profile_selected = "mm10"
                    params.pooledEI_snps_file = "/private/common/Software/AEI/RNAEditingIndex1.1/RNAEditingIndexer/Resources/SNPs/MusMusculus/ucscMM10CommonGenomicSNPs142.bed.gz"
                }
                standard {
                    params.profile_selected = ""
                    params.pooledEI_snps_file = ""

                }


            '''
            .stripIndent()

}

// if(!params.RID_or_bed6_file && !params.dynamic_regions) {
//   println "you must set RID_file or pool_bed6_file"
//   System.exit(1)
// }


process pool_regions_index {
  input:
    path cmpileups
    path pooling_result_dir
    path RID_or_bed6_file
    path snps_file
  output:
    path pooling_result_dir
  script:
    """
    # create RID file if wasnt given
    if [ $params.use_RID_file = true ]; then
      cat $RID_or_bed6_file > ./RID_file.csv
    else
      $params.which_createRID_command -i \$(readlink -f $RID_or_bed6_file) -o ./RID_file.csv
    fi
    # create group file if wasnt given
    if [ -z $params.group_file_csv ]; then
      echo "sample,group" > ./group_file.txt
      for d in \$(ls $cmpileups); do echo -e "\$d,\$d" >> ./group_file.txt; done
    else 
      cat $params.group_file_csv > ./group_file.txt
    fi
    # run pooling script
    $params.which_pooling_script_command -i $cmpileups -o \$(readlink -f $pooling_result_dir) --group_file ./group_file.txt -rid ./RID_file.csv --snps \$(readlink -f $snps_file) -s mpileup.cmpileup -p $params.pool_CMP_max_parallell
    """
}
  // $params.which_pooling_script_command -i $cmpileups -o \$(readlink -f $pooling_result_dir) --group_file ./group_file.txt -rid ./RID_file.csv -s ${params.cmp_regions_name}_mpileup.cmpileup -p $params.pool_CMP_max_parallell"""


// include { CMPILEUP_PIPELINE } from "${params.levanon_lab_pipelines_dir}/Subpipelines/cmpileup.nf"
include { CMPILEUP_Dynamic_regions_PIPELINE } from "${params.levanon_lab_pipelines_dir}/Subpipelines/cmpileup.nf" addParams(run_cmp_summarize : false)


// this mode dosent support RID file
workflow POOLED_EI_Dynamic_regions_PIPELINE {
    take:
    All_bams_dir_path
    regions
    bed6_file
    result_dir
    star_finished
    main:
    // CMPILEUP_PIPELINE will make a bam channle from the dir
    cmp_res = CMPILEUP_Dynamic_regions_PIPELINE(All_bams_dir_path, regions, result_dir, star_finished)
    pooled_index_dir = pool_regions_index(cmp_res.last(), result_dir, bed6_file, params.pooledEI_snps_file)
    emit:
    pooled_index_dir
}

// workflow POOLED_EDITING_INDEX_PIPELINE {
//     take:
//     All_bams_dir_path
//     result_dir
//     star_finished
//     main:
//     // CMPILEUP_PIPELINE will make a bam channle from the dir
//     cmp_res = CMPILEUP_PIPELINE(All_bams_dir_path, result_dir, star_finished)
//     pooled_index_dir = pool_regions_index(cmp_res.first(), params.RID_or_bed6_file,result_dir)
//     emit:
//     pooled_index_dir
// }





 

workflow {
    if(params.help){
        helpMessage()
    } else {
        POOLED_EI_Dynamic_regions_PIPELINE(params.bams_dir,params.regions,params.pooledEI_result_dir, true).view()
        // POOLED_EI_Dynamic_regions_PIPELINE(params.bams_dir,params.RID_or_bed6_file,params.pooledEI_result_dir, true).view()
    }
}