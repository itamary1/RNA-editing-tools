

/*

// file name without path and suffixes
regions_without_suffix=(new File(params.one_regions_file)).name
while(regions_without_suffix.contains('.')){
    regions_without_suffix = regions_without_suffix.lastIndexOf('.').with {it != -1 ? regions_without_suffix[0..<it] : regions_without_suffix}
}


// name of result dir 
params.both_strands_resDir="${params.Eindex_result_dir}/${regions_without_suffix}"


// in the below blocks we will save the regions and the refseq files as first and second according to params.plusMinusWords[0],params.plusMinusWords[1]
if(params.one_regions_file.contains(params.plusMinusWords[0])) {
    params.regions_first_file=params.one_regions_file
    params.regions_second_file=params.one_regions_file.replace(params.plusMinusWords[0],params.plusMinusWords[1])
} else {
    if(params.one_regions_file.contains(params.plusMinusWords[1])) {
    params.regions_second_file=params.one_regions_file
    params.regions_first_file=params.one_regions_file.replace(params.plusMinusWords[1],params.plusMinusWords[0])
    } else {
        if(!params.help){
        println "file one_regions_file should contain ${params.plusMinusWords[0]} or ${params.plusMinusWords[1]}, exiting.."
        exit 1
        }
    }
}

// same for refseq files
if(params.one_refseq_file.contains(params.plusMinusWords[0])) {
    params.refseq_first_file=params.one_refseq_file
    params.refseq_second_file=params.one_refseq_file.replace(params.plusMinusWords[0],params.plusMinusWords[1])
} else {
    if(params.one_refseq_file.contains(params.plusMinusWords[1])) {
    params.refseq_second_file=params.one_refseq_file
    params.refseq_first_file=params.one_refseq_file.replace(params.plusMinusWords[1],params.plusMinusWords[0])
    } else {
        if(!params.help){
        println "file one_refseq_file should contain ${params.plusMinusWords[0]} or ${params.plusMinusWords[1]}, exiting.."
        exit 1
        }
    }
}

//same for RID file

if(params.one_RID_file.contains(params.plusMinusWords[0])) {
    params.RID_first_file=params.one_RID_file
    params.RID_second_file=params.one_RID_file.replace(params.plusMinusWords[0],params.plusMinusWords[1])
} else {
    if(params.one_RID_file.contains(params.plusMinusWords[1])) {
    params.RID_second_file=params.one_RID_file
    params.RID_first_file=params.one_RID_file.replace(params.plusMinusWords[1],params.plusMinusWords[0])
    } else {
        if(!params.help){
        println "file one_RID_file should contain ${params.plusMinusWords[0]} or ${params.plusMinusWords[1]}, exiting.."
        exit 1
        }
    }
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
include { POOLED_EDITING_INDEX_PIPELINE as FIRST_EI} from "./editing_index_pooled_regions.nf" addParams(regions: "$params.regions_first_file",refseq_file: "$params.refseq_first_file", genome_regions_index : "${params.regions_first_file}.GenomeIndex.jsd", alu_index : false, RID_file : RID_first_file)
include { POOLED_EDITING_INDEX_PIPELINE as SECOND_EI} from "./editing_index_pooled_regions.nf" addParams(regions: "$params.regions_second_file",refseq_file: "$params.refseq_second_file", genome_regions_index : "${params.regions_second_file}.GenomeIndex.jsd", alu_index : false, RID_file : RID_second_file)
include { COMBINE_WAITING } from "${params.levanon_lab_pipelines_dir}/Subpipelines/combine_waiting.nf"

    
    
    
workflow POOLED_Eindex_plusMinus_PIPELINE {
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
        POOLED_Eindex_plusMinus_PIPELINE(params.indir, true)
        POOLED_Eindex_plusMinus_PIPELINE.out.view()
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
*/

// process split_to_plusMinus {
//   input:
//     path 'bed6file_regions.bed'
//     path 'refseq_file.bed'
//   output:
//     path './bed6file_regions.plus.bed' emit: plus_regions
//     path './bed6file_regions.minus.bed' emit: minus_regions
//   script:
//   """
//   $params.read_regions_command bed6file_regions.bed | awk 'BEGIN { FS = "\t"}; \$6=="+" {print \$0}' | bedtools sort | bedtools merge -c 4,5,6 -o distinct > ./bed6file_regions.plus.bed
//   $params.read_regions_command bed6file_regions.bed | awk 'BEGIN { FS = "\t"}; \$6=="-" {print \$0}' | bedtools sort | bedtools merge -c 4,5,6 -o distinct > ./bed6file_regions.minus.bed
//   """
// }