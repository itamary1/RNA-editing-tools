
nextflow.enable.dsl=2
params.help=false
params.have_base_config=false
if((!params.have_base_config) && (!params.help)){
 println "please attach a base config file"
 System.exit(1) 
}



def helpMessage() {
    log.info '''\

            note it will compress fastq at the end unless U set the flag to false

            run editing analysis from sra-acc list  
            ===================================
            nohup nextflow -c /home/alu/twerski/Scripts/Nextflow/Levanon_lab_NEXTFLOW_PIPELINE/Configs/Docker/full_levanon_pipeline.nf.docker.config -bg run /home/alu/twerski/Scripts/Nextflow/Levanon_lab_NEXTFLOW_PIPELINE/full_levanon_pipeline.nf -profile hg38 --use_existing_fastq false --ACC_list /home/alu/twerski/Scripts/Nextflow/training_data/small_srr.RL75.Human.txt --genome_length 75 --project_dir $PWD &> run.out.txt &




            >> important params
            directory for your pipeline files and results, if wont set it wil be $PWD
            params.project_dir="$launchDir"
            params.genome_length=''
            // if download - set acc list
            params.use_existing_fastq=false
            params.ACC_list=''
            // if not download - set fastqs input dir
            // directory containing fastq files >> use if use_existing_fastq was set to true <<
            params.indir=''



            //optional prameters depend to combine with indir:
            params.fastq_suffix='fastq'
            params.mates_patt="_{1,2}"
            params.fastq_pat="${params.indir}/*${params.mates_patt}.${params.fastq_suffix}"

            >>>>>general params
            params.raw_data_dir="${params.project_dir}/Raw_data"
            params.singleEnd=false
            params.results_dir="${params.project_dir}/Results"

            >>>>optional steps
            params.lzma_original_fastq=true
            if to add step of trimming and run FASTQC again
            note that set to this flag to true will result with 3 fastq dirs and matching FASTQC - orig/dedup/dedup_trim
            params.add_trimStep=false
            // if you want to trim in addtion to dedups of fastp (in one step) - note to set num ouf bp to trim below
            params.trim_in_dedup=false

            choose how many bases to trim from fromt and tail
            params.trimstep_basesF=10
            params.trimstep_basesT=10

            
            params.annotation='hg38'
            >> salmon params
            you need to override in case of not hg38 and if you want to run salmon:
            params.tx2id_geneMap = '/private/dropbox/Salmon_1.4.0/salmon_index/hg38/gencode_v32.transcriptToGeneID.tab'
            params.transcripts_index == '/private/dropbox/Salmon_1.4.0/salmon_index/hg38'
            
            >> STAR is supporting hg38 and mm10 with several lengthes, in case of other genome inedexes override
            params.star_genome
            
            >>>>flags and params
            params.singleEnd=false
            params.run_salmon=true
            params.run_AEI=true
            params.compress_fastqs=true

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
            you can also set othe specific params rlevats to the subpipelines, (TODO write the list in that message)

            '''
            .stripIndent()

}



include { DOWNLOAD_AND_PREPROCESS_PIPELINE } from './downloadAndPreprocess.nf' addParams(project_dir: "$params.project_dir")
include { ANALIZE_EDITING_PIPELINE } from './analizeEditing_on_cleanData.nf' addParams(project_dir: "$params.project_dir")


// not sure/remember why I need that, but its work
// process FAKE_PROCCESS {
//     input:
//     val(all_finished)
//     script:
//     finished = true
//     """
//     echo finished > ${params.project_dir}/finished.txt
//     """
// }

workflow FULL_LEVANON_PIPELINE {
    take:
    fastq_source
    main:
    cleaned_fastq = DOWNLOAD_AND_PREPROCESS_PIPELINE(fastq_source)
    ANALIZE_EDITING_PIPELINE(cleaned_fastq)
    emit:
    bams_dirs_ch = ANALIZE_EDITING_PIPELINE.out.bams_dirs_ch
    pipeline_finished = ANALIZE_EDITING_PIPELINE.out.everything_finished
}


workflow {
    if(params.help){
        helpMessage()
    } else {
        if (params.use_existing_fastq)
            FULL_LEVANON_PIPELINE(params.fastq_pat)
        else
            FULL_LEVANON_PIPELINE(params.ACC_list)
    }
}

workflow.onComplete {
    log.info("""
    Complete: Workflow Ended
    """)
}

workflow.onError {
    log.info("""
    Error: Workflow Ended With Error
    """)
}