nextflow.enable.dsl=2
params.help=false
params.have_base_config=false
if((!params.have_base_config) && (!params.help)){
 println "please attach a base config file"
 System.exit(1) 
}



def helpMessage() {
    log.info '''\

    
            compress fastqs workflow
            ===================================                   
            set yout indir - params.indir"
            !!!!docker supporting gzip only!!!!!

            >>> other params
            params.fastq_suffix='fastq'
            params.mates_patt='_{1,2}'
            params.fastq_pat="${params.indir}/*${params.mates_patt}.${params.fastq_suffix}"
            params.singleEnd=false
            '''
            .stripIndent()

}




process compress_files {
    tag "compress $acc"
    input:
        tuple val(acc), path(fastqs)
        val all_finished
    output:
        val fin
    script:
        fin = true
        def single = fastqs instanceof Path
        def read1 = !single ? /"${fastqs[0]}"/ : /"${fastqs}"/
        def read2 = !single ? /"${fastqs[1]}"/ : ''
        """
        $params.compress_command \$(readlink -f ${read1})
        if [ $single = false ]; then
            $params.compress_command \$(readlink -f ${read2})
        fi
        """
}


workflow COMPRESS_FASTQs {
    take:
    fastq_ch
    finished
    main:
    all_finished=compress_files(fastq_ch,finished)
    emit:
    all_finished
}

workflow{
    if(params.help){
        helpMessage()
    } else {
    fastq_ch=Channel.fromFilePairs( params.fastq_pat, size: params.singleEnd ? 1 : 2 )
        .ifEmpty { System.exit(1), "Cannot find any reads matching: ${params.fastq_pat}\nIf this is single-end data, please specify --singleEnd on the command line." }
    COMPRESS_FASTQs(fastq_ch,[true, true])
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