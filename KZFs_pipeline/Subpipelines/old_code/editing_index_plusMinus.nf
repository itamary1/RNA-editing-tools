nextflow.enable.dsl=2

/*******************
 config example at " /private8/Projects/itamar/ZNF/jose/Differentiated_PA1/Results/EditingIndex_perRegionKeepCMP/KRABS_locations_gencode/editing_index.config"
******************/


params.help=false


// if((!(params.one_regions_file && params.one_regions_file)) && !params.help){
//     println "please set one regions_file and one refseq_file"
//     exit 1
// }


// def helpMessage() {
//     log.info '''\

//             this is version of runing editing index for case of running in the both strands:
//             As every region in regular editing index run will be tested on one strand according to refseq annotation.
//             Here we allow you to seprate the editing index run to two different runs according to the strand you choose.
//             you need for the 2 regions files and 2 pseudo refseq files
//             ===================================
//             set:
//             params.indir=''
//             params.one_regions_file=""
//             params.one_refseq_file=""
//             // you need to have two regions and two refseq files.
//             //  the name of each file should contain word from params.plusMinusWords describe the strand of the file 
//             // the pipeline will access the opoosite strand file by replacing the word to the second option 
//             // for example: you give the file "KRABS_locations_gencode.bed3.minus.bed"

//             output will be writen to ${params.proj_resDir}/${regions_without_suffix} 
//             - regions_without_suffix is your regions file name
            

//             plusMinusWords=["plus", "minus"]

//             Eindex_result_dir="${launchDir}/Editing_index/"
//             indir=''
//             // !!! remeber to set your profile for othe params for editing_index.nf.docker.config
//             // index run params: 
//             PMkeep_cmpileup=false
//             PMgenome_file=''
//             PMbams_suffix='.Aligned.sortedByCoord.out.bam'
//             PMper_region_output=false
//             PMper_sample_output=false
//             PMsnp_file=''
//             PMexpression_file=''
//             // index will run only on path that contain that string -
//             PMmust_contain=''

//             '''
//             .stripIndent()

// }







// // file name without path and suffixes
// regions_without_suffix=(new File(params.one_regions_file)).name
// while(regions_without_suffix.contains('.')){
//     regions_without_suffix = regions_without_suffix.lastIndexOf('.').with {it != -1 ? regions_without_suffix[0..<it] : regions_without_suffix}
// }


// // name of result dir 
// params.both_strands_resDir="${params.Eindex_result_dir}/${regions_without_suffix}"


// // in the below blocks we will save the regions and the refseq files as first and second according to params.plusMinusWords[0],params.plusMinusWords[1]
// if(params.one_regions_file.contains(params.plusMinusWords[0])) {
//     params.regions_first_file=params.one_regions_file
//     params.regions_second_file=params.one_regions_file.replace(params.plusMinusWords[0],params.plusMinusWords[1])
// } else {
//     if(params.one_regions_file.contains(params.plusMinusWords[1])) {
//     params.regions_second_file=params.one_regions_file
//     params.regions_first_file=params.one_regions_file.replace(params.plusMinusWords[1],params.plusMinusWords[0])
//     } else {
//         if(!params.help){
//         println "file one_regions_file should contain ${params.plusMinusWords[0]} or ${params.plusMinusWords[1]}, exiting.."
//         exit 1
//         }
//     }
// }

// // same for refseq files
// if(params.one_refseq_file.contains(params.plusMinusWords[0])) {
//     params.refseq_first_file=params.one_refseq_file
//     params.refseq_second_file=params.one_refseq_file.replace(params.plusMinusWords[0],params.plusMinusWords[1])
// } else {
//     if(params.one_refseq_file.contains(params.plusMinusWords[1])) {
//     params.refseq_second_file=params.one_refseq_file
//     params.refseq_first_file=params.one_refseq_file.replace(params.plusMinusWords[1],params.plusMinusWords[0])
//     } else {
//         if(!params.help){
//         println "file one_refseq_file should contain ${params.plusMinusWords[0]} or ${params.plusMinusWords[1]}, exiting.."
//         exit 1
//         }
//     }
// }

process split_to_plusMinus {
  input:
    path 'bed6file_regions.bed'
    path 'refseq_file.bed'
  output:
    './bedfile_regions.plus.bed' emit: plus_regions
    './bedfile_regions.minus.bed' emit: minus_regions
  script:
  """
  $params.read_regions_command bed6file_regions.bed | awk 'BEGIN { FS = "\t"}; \$6=="+" {print \$0}' | bedtools sort | bedtools merge > ./bedfile_regions.plus.bed
  $params.read_regions_command bed6file_regions.bed | awk 'BEGIN { FS = "\t"}; \$6=="-" {print \$0}' | bedtools sort | bedtools merge > ./bedfile_regions.minus.bed
  $params.read_refseq_command refseq_file.bed |  awk 'BEGIN { FS = "\t"}; \$6=="+" {print \$0}' > ./refseq_file.plus.bed
  $params.read_refseq_command refseq_file.bed |  awk 'BEGIN { FS = "\t"}; \$6=="-" {print \$0}' > ./refseq_file.minus.bed
  """
}

process create_dirs {
    input:
        path EI_dir
        val both_strand_dir
    output:
        env f_outdir_out
        env s_outdir_out
    script:
    """
    echo 
    f_outdir=${both_strand_dir}/plus
    s_outdir=${both_strand_dir}/minus
    mkdir -p \$f_outdir
    mkdir -p \$s_outdir
    f_outdir_out=\$f_outdir
    s_outdir_out=\$s_outdir
    """
}

// EI_PIPELINE take: All_bams_dir_path,result_dir,star_finished emits:  RUN_EINDEX.out
include { EI_PIPELINE as FIRST_EI} from "${params.levanon_lab_pipelines_dir}/Subpipelines/editing_index.nf" addParams(regions: "$params.regions_first_file",refseq_file: "$params.refseq_first_file", genome_regions_index : "${params.regions_first_file}.GenomeIndex.jsd", alu_index : false)
include { EI_PIPELINE as SECOND_EI} from "${params.levanon_lab_pipelines_dir}/Subpipelines/editing_index.nf" addParams(regions: "$params.regions_second_file",refseq_file: "$params.refseq_second_file", genome_regions_index : "${params.regions_second_file}.GenomeIndex.jsd", alu_index : false)
include { COMBINE_WAITING } from "${params.levanon_lab_pipelines_dir}/Subpipelines/combine_waiting.nf"

    
    
    
workflow Eindex_plusMinus_PIPELINE {
    take:
    bams_dir
    readyToGo
    main:
    create_dirs(params.Eindex_result_dir,params.both_strands_resDir)
    first_outdir=create_dirs.out[0]
    second_outdir=create_dirs.out[1]
    F_finished = FIRST_EI(bams_dir,first_outdir, readyToGo)
    S_finished = SECOND_EI(bams_dir,second_outdir, readyToGo)
    all_finished = COMBINE_WAITING(F_finished,S_finished)
    emit:
    all_finished
}



workflow {
    if (params.help)
        helpMessage()
    else
        Eindex_plusMinus_PIPELINE(params.indir, true)
        Eindex_plusMinus_PIPELINE.out.view()
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
