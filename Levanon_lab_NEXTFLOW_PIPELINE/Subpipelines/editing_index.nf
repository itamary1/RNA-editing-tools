nextflow.enable.dsl=2

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
            RUN RNA EDITING INDEXER
            ===================================

            commnad example:
            config=/home/alu/twerski/Scripts/Nextflow/Levanon_lab_NEXTFLOW_PIPELINE/Configs/Docker/SubP_configs/editing_index.nf.docker.config
            script=/home/alu/twerski/Scripts/Nextflow/Levanon_lab_NEXTFLOW_PIPELINE/Subpipelines/editing_index.nf
            bams_dir=/private10/Projects/Itamar/check_Lab_pipline/try_FUP/Raw_data/STAR
            nohup nextflow -bg -c $config run $script -profile hg38 --bams_dir $bams_dir --Eindex_result_dir $PWD/Results/AEI &> run_EI.out.txt &

            by defualt it will run on alu index, if you want to change it -
            set regions, alu_index(=false)
            and probably set:
            Eindex_result_dir[="${launchDir}/AEI/"]

            set -profile hg38/mm10
            otherwise you should set:
            params.genome ('UserProvided')
            params.genome_file 
            params.refseq_file 
            params.expression_file
            params.snp_file 


            run in docker: lite docker containing no resources(genome, expression file ..) in which you must specify all the resources.
            I already create profile for local hg38 and mm10 resources - run with "-profile "
            if somone have time to add more profiles...
            there is a full version of docker which support built in files for several popular assemblies, but there is no support for it yet(becaouse docker mounting constraint, and I dont think its important/usfull).

            all params and their defaults:
            bams_dir =''
            Eindex_result_dir="${launchDir}/AEI/"
            //run in lite version - use local resources files instead of inner docker's
            // currently not supporting run on full docker (docker_lite=false)
            docker_lite=true
            stranded=false
            singleEnd=false
            keep_cmpileup=false
            genome_file=''
            bams_suffix='.Aligned.sortedByCoord.out.bam'
            per_region_output=false
            per_sample_output=false
            regions=''
            snp_file=''
            refseq_file=''
            expression_file=''
            // index will run only on path that contain that string -
            must_contain=''
            which_RNAEditingIndex='RNAEditingIndex'
            additional_params = ''
            alu_index=true
            // if U want to run on many regions files - set to true and set regions to be a pattern like '/path/to/files/*.bed'
            // if U run regions='/path/to/files/*.bed' ; nextflow -c .. run .. --regions "$regions" >>> remember to use quotes 
            run_on_multiRegionsFiles = false

            '''
            .stripIndent()

}

include { COMBINE_WAITING as COMB } from "${params.levanon_lab_pipelines_dir}/Subpipelines/combine_waiting.nf"

process build_EI_docker {
    maxForks 1
    input:
        val go
        path EI_Docker_file
    output:
        val finished
    script:
        finished = true
        """
            # if docker image not exist in the system - build it
            if [ -z "\$(docker images -q $params.editing_index_docker_image 2> /dev/null)" ]; then
                git clone https://github.com/a2iEditing/RNAEditingIndexer.git
                rm ./RNAEditingIndexer/Dockerfile
                cp $params.EI_Docker_file ./RNAEditingIndexer/Dockerfile
                cd ./RNAEditingIndexer/
                docker image build -t $params.editing_index_docker_image .
                cd ../ && rm -r -f ./RNAEditingIndexer/
            else
            echo "docker exist, skipping"; fi
        """
}

// some pathes needed for run in docker
// if(params.docker_run){
//     result_dirBaseName=(new File(params.Eindex_result_dir)).name
//     result_dir="/data/output/$result_dirBaseName"
// }



process RUN_EINDEX {
    tag "editing index $All_bams_dir_path_full"
    maxForks params.max_fork_RNAEditingIndex
    input:
        val All_bams_dir_path_full
        path All_bams_dir_path
        path result_dir
        val full_res_dir
        path expression_file
        path refseq_file
        path snp_file
        path regions
        val regions_full
        path genome_file
        val star_finished
        output:
        path './index_result'
    script:
        // take care of all the flags
        
        if(params.stranded){
            stranded_flag='--stranded'
            if(params.singleEnd)
                paired_flag=''
            else
                paired_flag='--paired_end'
        }
        else {
            stranded_flag=''
            paired_flag=''
        }

        // most of the flags are mandatory in the lite docker, but I did it optional to support the full docker - (which is maybe not fully supported in other parts of the code)
        mustContain = params.must_contain ? "-s $params.must_contain" : ''
        per_region_flag = params.per_region_output ? '--per_region_output' : ''
        per_sample_flag = params.per_sample_output ? '--per_sample_output' : ''
        keep_cmpileup_flag = params.keep_cmpileup ? '--keep_cmpileup' : ''
        genome_file_command = params.genome_file ? "-gf $params.genome_file" : ''
        regions_command = regions_full != params.EmptyPath ? "-rb $regions_full" : ''
        snp_command = (params.snp_file && params.snp_file != params.EmptyPath) ? "--snps $params.snp_file" : ''
        refseq_command = refseq_file != params.EmptyPath ? "--refseq $refseq_file" : ''
        //expression file is optional - could be empty path
        expression_command = (params.expression_file && params.expression_file != params.EmptyPath) ? "--genes_expression $params.expression_file" : ''
        additional_command = params.editingIndex_additional_params ? " -a $params.editingIndex_additional_params " : ''

        """
        regions_command="${regions_command}"
        mkdir -p ${full_res_dir}/logs
        mkdir -p ${full_res_dir}/cmpileups
        mkdir -p ${full_res_dir}/summary
        if [ "$regions_command" ]; then
            num_cols=\$(head -1 $regions_full | awk '{print NF}')
            # if regions need fixing  - fix it
            if [ $params.mergeBED_forEI = true ] || [ \$num_cols != 3 ]; then
                cat $regions_full | bedtools sort | bedtools merge > ./EI_regions.merged.bed3.bed
                # update regions command
                regions_command="-rb \${PWD}/EI_regions.merged.bed3.bed"
            fi
        fi
        $params.which_RNAEditingIndex -d $All_bams_dir_path_full -f $params.bams_suffix $stranded_flag $mustContain $paired_flag $per_region_flag $per_sample_flag $keep_cmpileup_flag  -l ${full_res_dir}/logs -o ${full_res_dir}/cmpileups -os ${full_res_dir}/summary --follow_links --genome $params.genome $genome_file_command $refseq_command $snp_command ${expression_command}${additional_command} \$regions_command
        ln -s ${full_res_dir} ./index_result
        """
}

// to handle run EI on list of regions (list of files)-
// we accept only one input of regions (not two of path and val like the regular) and this input can be list/channel
// regions command was adapted
// output file will be inner directory for each region to handle collusions
process RUN_EINDEX_multi_regions {
    tag "editing index $All_bams_dir_path_full"
    maxForks params.max_fork_RNAEditingIndex
    input:
        val All_bams_dir_path_full
        path All_bams_dir_path
        path result_dir
        val full_res_dir
        path expression_file
        path refseq_file
        path snp_file
        each path(regions)
        path genome_file
        val star_finished
        output:
        path './index_result'
    script:
        // take care of all the flags
        
        if(params.stranded){
            stranded_flag='--stranded'
            if(params.singleEnd)
                paired_flag=''
            else
                paired_flag='--paired_end'
        }
        else {
            stranded_flag=''
            paired_flag=''
        }

        // most of the flags are mandatory in the lite docker, but I did it optional to support the full docker - (which is maybe not fully supported in other parts of the code)
        mustContain = params.must_contain ? "-s $params.must_contain" : ''
        per_region_flag = params.per_region_output ? '--per_region_output' : ''
        per_sample_flag = params.per_sample_output ? '--per_sample_output' : ''
        keep_cmpileup_flag = params.keep_cmpileup ? '--keep_cmpileup' : ''
        genome_file_command = params.genome_file ? "-gf $params.genome_file" : ''
        // regions_command = "-rb $regions_full" : ''
        snp_command = (params.snp_file && params.snp_file != params.EmptyPath) ? "--snps $params.snp_file" : ''
        refseq_command = refseq_file != params.EmptyPath ? "--refseq $refseq_file" : ''
        //expression file is optional - could be empty path
        expression_command = (params.expression_file && params.expression_file != params.EmptyPath) ? "--genes_expression $params.expression_file" : ''
        additional_command = params.editingIndex_additional_params ? " -a $params.editingIndex_additional_params " : ''
        // for being able to run on list of regions in one result dir without collusions 
        // we create and use subdir for each region as full_res_dir
        as_str="$regions"
        no_suff=as_str.indexOf('.bed').with {it != -1 ? as_str[0..<it] : as_str}
        full_res_dir="${full_res_dir}/${no_suff}"
        """
        mkdir -p ${full_res_dir}/logs
        mkdir -p ${full_res_dir}/cmpileups
        mkdir -p ${full_res_dir}/summary
        $params.which_RNAEditingIndex -d $All_bams_dir_path_full -f $params.bams_suffix $stranded_flag $mustContain $paired_flag $per_region_flag $per_sample_flag $keep_cmpileup_flag  -l ${full_res_dir}/logs -o ${full_res_dir}/cmpileups -os ${full_res_dir}/summary --follow_links --genome $params.genome $genome_file_command $refseq_command $snp_command ${expression_command}${additional_command} -rb \$(readlink -f $regions)
        ln -s ${full_res_dir} ./index_result
        """
}

// workflow EI_PIPELINE_dynamic_regionsAndRefseq {
//     take:
//     All_bams_dir_path
//     regions
//     refseq_file
//     result_dir
//     star_finished
//     main:
//     bams_dir_ch=Channel.fromPath(All_bams_dir_path, type:'dir').ifEmpty { System.exit(1), "Cannot find your bams folder"}
//     RUN_EINDEX(All_bams_dir_path,\
//     bams_dir_ch,\
//     result_dir,\
//     params.genome_regions_index,\
//     params.expression_file,\
//     refseq_file\
//     params.snp_file,\
//     regions,\
//     params.genome_file,\
//     star_finished)
//     emit:
//     RUN_EINDEX.out
// }


workflow EI_PIPELINE_dynamicRegRef {
    take:
    All_bams_dir_path
    refseq_f
    regions_bed
    result_dir
    star_finished
    main:
    if((params.docker_run) && !(params.use_editing_index_hub_image)){
        docker_built = build_EI_docker(star_finished,params.EI_Docker_file)
        readyToGo=COMB(star_finished,docker_built)
    } else {
        readyToGo=star_finished
    }
    bams_dir_ch=Channel.fromPath(All_bams_dir_path, type:'dir').ifEmpty { System.exit(1), "Cannot find your bams folder"}
    if (params.run_on_multiRegionsFiles){
        RUN_EINDEX_multi_regions(All_bams_dir_path,\
        bams_dir_ch,\
        result_dir,\
        result_dir,\
        params.expression_file,\
        refseq_f,\
        params.snp_file,\
        regions_bed,\
        params.genome_file,\
        readyToGo)
        res = RUN_EINDEX_multi_regions.out
    } else {
        RUN_EINDEX(All_bams_dir_path,\
        bams_dir_ch,\
        result_dir,\
        result_dir,\
        params.expression_file,\
        refseq_f,\
        params.snp_file,\
        regions_bed,\
        regions_bed,\
        params.genome_file,\
        readyToGo)
        res = RUN_EINDEX.out
    }
    emit:
    res
}

workflow EI_PIPELINE {
    take:
    All_bams_dir_path
    result_dir
    star_finished
    main:
    EI_PIPELINE_dynamicRegRef(All_bams_dir_path, params.refseq_file,params.regions,result_dir,star_finished)
    emit:
    EI_PIPELINE_dynamicRegRef.out
}

workflow {
    if(params.help){
        helpMessage()
    } else {
        if (params.run_on_multiRegionsFiles){
            if (params.alu_index) {
                println "please set alu_index to false or run_on_multiRegionsFiles to false"
                System.exit(1)
            }
            println "$params.regions"
            regions_ch = Channel.fromPath("$params.regions").ifEmpty { System.exit(1), "Cannot find bed files at $params.regions"}
            EI_PIPELINE_dynamicRegRef(params.bams_dir, params.refseq_file,regions_ch,params.Eindex_result_dir,true).view()
        } else {
            EI_PIPELINE(params.bams_dir,params.Eindex_result_dir, true).view()
        }
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