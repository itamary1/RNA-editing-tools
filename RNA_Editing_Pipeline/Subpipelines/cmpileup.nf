
params.help=false
params.have_base_config=false
if((!params.have_base_config) && (!params.help)){
 println "please attach a base config file"
 System.exit(1)
}


// if((!params.profile_selected) && (!params.help)){
// println "please select a profile with -profile <proile_name>"
//  System.exit(1) 
// }

// if((params.docker_run) && (!params.bams_dir) && (!params.help)){
//     println "you must set bams dir when runing with docker. exiting.."
//     System.exit(1)
// }

def helpMessage() {
    log.info '''\
            RUN CMpileup
            ===================================
            command example:
            bams_dir="/private10/Projects/Itamar/check_Lab_pipline/try_FUP/Raw_data/STAR"
            nohup nextflow -bg -c /home/alu/twerski/Scripts/Nextflow/Levanon_lab_NEXTFLOW_PIPELINE/Configs/Docker/SubP_configs/cmpileup.nf.docker.config run /home/alu/twerski/Scripts/Nextflow/Levanon_lab_NEXTFLOW_PIPELINE/Subpipelines/cmpileup.nf -profile hg38 --bams_dir $bams_dir --cmp_result_dir $PWD/Results/CMP &> run_CMP.out.txt &

            note that if you giving list/channel of regions, you shouldnt set params.cmp_regions_name because of names conflict

            params and defaults:
            stranded = false
            singleEnd = false
            bams_dir = ''
            cmp_result_dir = "${launchDir}/CMPileup"
            cmp_regions_name = ''
            which_CMP_python = 'python'
            which_CMP_PipelineManger = '/home/biodocker/GGPS/Session/PipelineManger.py'
            cmp_truncation = '0,0'
            cmp_additional_params = ''
            regions_name = 'unset'
            bams_suffix='.Aligned.sortedByCoord.out.bam'
            run_soft_intersect = false 
            conf_prefix= params.run_soft_intersect ? "${params.levanon_lab_pipelines_dir}/Resources/CMPileup/Docker/PileupKnown.SoftIntersect" : "${params.levanon_lab_pipelines_dir}/Resources/CMPileup/Docker/PileupKnown"
            cmp_conf_file = params.stranded ? "${params.conf_prefix}.stranded.conf" : "${params.conf_prefix}.stranded.conf"
            // can be bed 3 or more
            cmp_regions_bed = "${params.levanon_lab_pipelines_dir}/Resources/CMPileup/orshai_AG_only_published.bed6.bed"
            run_cmp_summarize = true
            which_cmp_summarize = "Rscript ${params.levanon_lab_pipelines_dir}/Resources/process_cmpileup.R"
            
            '''
            .stripIndent()

}


include { COMBINE_WAITING as COMB } from "${params.levanon_lab_pipelines_dir}/Subpipelines/combine_waiting.nf"

process build_CMP_docker {
  maxForks 1  
  input:
    val go
    path cmpileup_docker_image
  output:
    val finished
  script:
    finished = true
    """
      # if docker image not exist in the system - build it
      if [ -z "\$(docker images -q $params.cmpileup_docker_image 2> /dev/null)" ]; then
        docker image build -t $params.cmpileup_docker_image $params.cmpileup_Docker_file_dir
      else
        echo "docker exist, skipping"; fi
    """
}


// process RUN_CMP {
//     tag "cmpileup on $All_bams_dir_path_full"
//     // maxForks 1
//     input:
//         val All_bams_dir_path_full
//         path All_bams_dir_path
//         path result_dir
//         val full_res_dir
//         path cmp_conf_file
//         path genome_fa
//         path regions_bed
//         val star_finished
//     output:
//         path './cmp_result/cmpileups'
    
//     script:
//         """
//         awk '{OFS="\t"; print \$1,\$2,\$3}' $regions_bed > ${regions_bed}.bed3.bed
//         # create dirs for results and logs
//         mkdir -p ${full_res_dir}/logs
//         mkdir -p ${full_res_dir}/cmpileups
//         # save my currnet working directory
//         Mywd=\$(pwd)
//         $params.which_CMP_python $params.which_CMP_PipelineManger -t $params.cmp_truncation -c $cmp_conf_file -d $All_bams_dir_path_full -f $params.bams_suffix -o ${full_res_dir}/cmpileups -l ${full_res_dir}/logs --follow_links -a regions_coordinates_bed=\\'\${Mywd}/${regions_bed}.bed3.bed\\' regions_name=\\'${params.cmp_regions_name}\\' genome_fasta=\\'\$(readlink -f ${genome_fa})\\' bam_file_suffix=\\'${params.bams_suffix}\\' $params.cmp_additional_params
//         ln -s ${full_res_dir} ./cmp_result
//         """
// }

// in addition to each path(regions_bed) -
// in this process the logs are in local directory 
// if params.cmp_regions_name unset regions_name is same as regions_bed var removing bed* suffix
process RUN_CMP {
    tag "cmpileup on $All_bams_dir_path_full"
    maxForks 1
    input:
        val All_bams_dir_path_full
        path All_bams_dir_path
        path result_dir
        val full_res_dir
        path cmp_conf_file
        path genome_fa
        // U can run on many regions in one run
        each path(regions_bed)
        val star_finished
    output:
        path './cmp_result/cmpileups'
    
    script:
          paried_flag = (!(params.singleEnd) && params.stranded) ? "\\' is_paired_end=True\\'" : ''
        """
        awk '{OFS="\t"; print \$1,\$2,\$3}' $regions_bed > ${regions_bed}.bed3.bed
        # create dirs for results and logs
        mkdir ./logs
        mkdir -p ${full_res_dir}/cmpileups
        # save my currnet working directory
        Mywd=\$(pwd)
        # set name for your cmpileup files
        #if its empty - same as regions but without suffix
        if [ -z $params.cmp_regions_name ]; then
          regions_bed_no_suff=$regions_bed
          regions_bed_no_suff=\${regions_bed_no_suff%%.bed*}
        else
          regions_bed_no_suff=$params.cmp_regions_name
        fi
        # run cmpileup
        $params.which_CMP_python $params.which_CMP_PipelineManger -t $params.cmp_truncation -c $cmp_conf_file -d $All_bams_dir_path_full -f $params.bams_suffix -o ${full_res_dir}/cmpileups -l \${Mywd}/logs --follow_links -a regions_coordinates_bed=\\'\${Mywd}/${regions_bed}.bed3.bed\\' regions_name=\\'\${regions_bed_no_suff}\\' genome_fasta=\\'\$(readlink -f ${genome_fa})\\' bam_file_suffix=\\'${params.bams_suffix}\\'${paried_flag} ${params.cmp_additional_params}
        ln -s ${full_res_dir} ./cmp_result
        """
} 
    
process CMP_SUMMARY {
  input:
    path result_dir
    path cmp_finished
  output:
    stdout
  script:
  """
  $params.which_cmp_summarize -i $result_dir -s ${params.cmp_regions_name}_mpileup.cmpileup -p 20
  """
}

workflow CMPILEUP_Dynamic_regions_PIPELINE {
    take:
    All_bams_dir_path
    regions_bed
    result_dir
    star_finished
    main:
    if((params.docker_run) && !(params.use_cmpileup_hub_image)){
      dockerBuilt = build_CMP_docker(star_finished,params.cmpileup_Docker_file_dir)
      readyToGo = COMB(star_finished,dockerBuilt)
    } else {
      readyToGo=star_finished
    }
    bams_dir_ch=Channel.fromPath(All_bams_dir_path, type:'dir').ifEmpty { System.exit(1), "CMpileup: Cannot find your bams folder"}
    RUN_CMP(All_bams_dir_path,\
    bams_dir_ch,\
    result_dir,\
    result_dir,\
    params.cmp_conf_file,\
    params.genome_file,\
    regions_bed,\
    readyToGo)
    if(params.run_cmp_summarize)
        CMP_SUMMARY(result_dir,RUN_CMP.out)
        // params.cmp_snp_file,\
    emit:
    RUN_CMP.out
}
 
workflow CMPILEUP_PIPELINE {
    take:
    All_bams_dir_path
    result_dir
    star_finished
    main:
    CMPILEUP_Dynamic_regions_PIPELINE(All_bams_dir_path,params.cmp_regions_bed,result_dir,star_finished)
    emit:
    CMPILEUP_Dynamic_regions_PIPELINE.out
}

workflow {
    if(params.help){
        helpMessage()
    } else {
        CMPILEUP_PIPELINE(params.bams_dir,params.cmp_result_dir, true).view()
    }
}

workflow.onComplete {
    log.info("""
    Complete: Workflow Ended
    """)
}

workflow.onError {
    log.info("""
    Error: Workflow Ended With Error ${workflow.errorMessage}
    """)
}