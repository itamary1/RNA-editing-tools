
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




def helpMessage() {
    log.info '''\
            GET RNA EDITING INDEX ON POOLED REGIONS
            ==================================
            example command:
            OD=${PD}Scripts/Nextflow/training_data/trys/trypool && cd $OD 
            Conf=${PD}Scripts/Nextflow/Special_pipelines/Configs/Dockers/SubP_configs/editing_index_pooled_regions.nf.docker.config
            Script=${PD}Scripts/Nextflow/Special_pipelines/Subpipelines/editing_index_pooled_regions.nf
            Reg=${PD}Scripts/Nextflow/Special_pipelines/Resources/KZFs/CDS/inverted_KRABS_CDS.gencode.bed6.bed
            nextflow -c $Conf run $Script --bams_dir ${PD}Scripts/Nextflow/training_data/example_bam/ --pooledEIregions $Reg --pooledEI_result_dir $OD --pooledEIbed6_file $Reg  -profile hg38


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
                which_python_39 = 'python3'
                you shloud set profile -hg38/mm10
                or set the SNPs path:
                params.pooledEI_snps_file
                
            '''
            .stripIndent()

}




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




 

workflow {
    if(params.help){
        helpMessage()
    } else {
        POOLED_EI_Dynamic_regions_PIPELINE(params.bams_dir,params.pooledEIregions,params.pooledEIbed6_file,params.pooledEI_result_dir, true).view()
    }
}