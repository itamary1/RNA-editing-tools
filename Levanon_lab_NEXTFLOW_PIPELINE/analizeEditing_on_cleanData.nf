nextflow.enable.dsl=2
params.help=false
params.have_base_config=false
if((!params.have_base_config) && (!params.help)){
 println "please attach a base config file"
 System.exit(1) 
}

// if((!params.help) && (!params.docker_run) && params.run_known) {
//  println "runing known localy is not supported yet\n exiting.."
//  System.exit(1) 
// }


def helpMessage() {
    log.info '''\

    
            Analize RNA Editing on Clean Fastq PIPELINE => will run Salmon,Star,AluEditingIndex and compress the Fastq at the End
            ===================================     
            nohup nextflow -bg -c /home/alu/twerski/Scripts/Nextflow/Levanon_lab_NEXTFLOW_PIPELINE/Configs/Docker/analizeEditing_on_cleanData.nf.docker.config run /home/alu/twerski/Scripts/Nextflow/Levanon_lab_NEXTFLOW_PIPELINE/analizeEditing_on_cleanData.nf -profile hg38 --genome_length 75 --indir /private10/Projects/Itamar/check_Lab_pipline/tryDAP/Raw_data/Fastp_out_fastqs &> run.out.txt 



            directory for your pipeline results, if wont set it wil be $PWD
            params.project_dir="$launchDir"
            you must select read length at 
            params.genome_length=''
            
            >> salmon params
            you need to override in case of not hg38 and if you want to run salmon:
            params.tx2id_geneMap = '/private/dropbox/Salmon_1.4.0/salmon_index/hg38/gencode_v32.transcriptToGeneID.tab'
            params.transcripts_index == '/private/dropbox/Salmon_1.4.0/salmon_index/hg38'
            
            >> STAR is supporting hg38 and mm10 withj several lengthes, in case of other genome inedexes override
            params.star_genome
            
            >>>>flags and params
            params.singleEnd=false
            params.run_salmon=true
            params.run_AEI=true
            params.compress_fastqs=true
            //note! runing known is based on cmpileup docker, and is not supported localy
            params.run_known

            >> pathes
            params.fastq_suffix='fastq'
            params.mates_patt = "_{1,2}"
            params.fastq_pat="${params.fastqs_dir}/*${params.mates_patt}.${params.fastq_suffix}"
            params.results_dir="${params.project_dir}/Results"
            params.salmon_dir="${params.results_dir}/Salmon"
            params.star_dir="${params.project_dir}/Raw_data/STAR"
            // alu editing index results dir
            params.Eindex_result_dir="${params.results_dir}/AEI"

            // if set to false it will use lzma.
            // compressing command can be override completely by params.compress_command=
            params.useGzip=true
            '''
            .stripIndent()

}


process make_project_dirs {
  input:
    path project_dir
  output:
    val finished
  script:
  finished = true
  """
  mkdir -p ${params.results_dir}
  mkdir -p ${params.raw_data_dir}
  mkdir -p ${params.editing_index_dir}
  sleep 1
  """
}


// process add_acc {
//   input:
//     path bams_dir
//   output:
//     tuple val(acc), path(bams_dir)
//   script:
//   as_str = "$bams_dir"
//   acc = as_str.lastIndexOf('/').with {it != -1 ? as_str[it..<as_str.length()] : as_str}
//   """
//   """
// }


//you also should look in the salmon index properties!!!
// if((params.annotation!='hg38' && params.transcripts_index == '/private/dropbox/Salmon_1.4.0/salmon_index/hg38' && params.run_salmon) && !params.help ) {
//     println "please override salmon-index in params.transcripts_index and params.tx2id_geneMap"
//     System.exit(1)
// }



// take fastq_ch,outdir
include { SALMON_PIPELINE } from './Subpipelines/salmon.nf' 
if(!params.bams_as_input){
 // take: fastq_ch output_dir emit sorted_indexed_bams_ch
    include { STAR_PIPELINE } from './Subpipelines/star.nf'
 }
include { EI_PIPELINE as AEI } from './Subpipelines/editing_index.nf' // take  All_bams_dir_path result_dir star_finished   emit:RUN_EINDEX.out

include { EI_PIPELINE as UTR3EI } from './Subpipelines/editing_index.nf' addParams(regions : "${params.UTR3_regions}") // take  All_bams_dir_path result_dir star_finished   emit:RUN_EINDEX.out


include { CMPILEUP_PIPELINE as KNOWN } from './Subpipelines/cmpileup.nf' // take  All_bams_dir_path result_dir star_finished   emit:RUN_CMP.out

include { COMBINE_WAITING as COMB1 ; COMBINE_WAITING as COMB2; COMBINE_WAITING as COMB3; COMBINE_WAITING as COMB4; COMBINE_WAITING as COMB5; COMBINE_WAITING as COMBX } from './Subpipelines/combine_waiting.nf' //    take: waiting_one waiting_two emit: both_finished

include { COMPRESS_FASTQs } from './Subpipelines/compress_fastqs.nf' //    take: fastq_ch,finished



// TODO
// include { HE_PIPELINE  } from




workflow ANALIZE_EDITING_PIPELINE {
    take:
    fastq_or_bams_ch
    main:
    //create project dirs
    dir_created = make_project_dirs(params.project_dir)
    if(params.bams_as_input) {
        bams_dir=params.indir
        bams_dirs_ch = fastq_or_bams_ch
        all_bams_dirs_finished = dir_created
        salmon_finished = true
    } else {
        bams_dir=params.star_dir
        // run salmon
        if(params.run_salmon)
        salmon_finished = SALMON_PIPELINE(fastq_or_bams_ch,params.salmon_dir,dir_created).collect()
        else
        salmon_finished = true
        // align fastqs to sorted indexed bams
        bams_dirs_ch = STAR_PIPELINE(dir_created,fastq_or_bams_ch,params.star_dir)
        all_bams_dirs_finished = bams_dirs_ch.collect()      
    }
    // set var which marks that all proccess using fastq was finished
    fastq_finished = COMB1(salmon_finished,all_bams_dirs_finished)
    if(params.run_AEI){
        aei_res=AEI(bams_dir,params.AEI_dir,all_bams_dirs_finished)
    } else {
        aei_res =true
    }
    if(params.run_3UTR){
        // wait for AEI to prevent building docker by both
        ready_UTREI=COMBX(all_bams_dirs_finished,aei_res)
        UTR3_res=UTR3EI(bams_dir,params.EI3UTR_dir,ready_UTREI)
    } else {
        UTR3_res =true
    }
    everything_finished0=COMB2(aei_res,fastq_finished)
    everything_finished1=COMB3(UTR3_res,everything_finished0)
    if(params.run_known){
        cmp_res=KNOWN(bams_dir,params.known_dir,all_bams_dirs_finished)
    } else {
        cmp_res = true
    }
    everything_finished2=COMB4(cmp_res,everything_finished1)
    if(params.compress_fastqs && !params.bams_as_input){
        c_f=COMPRESS_FASTQs(fastq_or_bams_ch, fastq_finished)
    } else {
        c_f=true
    }
    everything_finished=COMB5(everything_finished2,c_f)
    emit:
    bams_dirs_ch
    everything_finished
}


workflow {
    if (params.help) {
        helpMessage()
    } else {
        if(params.bams_as_input) {
        println "bams as input"
        // create channel of tupels :[accession_number(directory name), path_object]
        bams_dirs_ch = Channel.fromPath(params.bams_full_patt, type: 'dir' ).map{ bam_d -> tuple("$bam_d".lastIndexOf('/').with {it != -1 ? "$bam_d"[it+1..<"$bam_d".length()] : "$bam_d"}, bam_d) }.ifEmpty { System.exit(1), "Cannot find any bams matching: ${params.bams_full_patt}\nplease make sure you set indir and bams_patt properly" }
        ANALIZE_EDITING_PIPELINE(bams_dirs_ch)
        ANALIZE_EDITING_PIPELINE.out.bams_out_dirs_ch.view()
        // ANALIZE_EDITING_PIPELINE.out.everything_finished.view()
        } else {
        println "fastq input"
        cleaned_fastq_ch =  Channel.fromFilePairs(params.fastq_pat, size: params.singleEnd ? 1 : 2 ).ifEmpty { System.exit(1), "Cannot find any reads matching: ${params.fastq_pat}\nIf this is single-end data, please specify --singleEnd on the command line." }
        ANALIZE_EDITING_PIPELINE(cleaned_fastq_ch)
        ANALIZE_EDITING_PIPELINE.out.bams_dirs_ch.view()
        // ANALIZE_EDITING_PIPELINE.out.everything_finished.view()
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