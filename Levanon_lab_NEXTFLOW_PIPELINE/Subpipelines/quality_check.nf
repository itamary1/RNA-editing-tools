nextflow.enable.dsl=2
params.help=false
params.have_base_config=false
if((!params.have_base_config) && (!params.help)){
 println "please attach a base config file"
 System.exit(1) 
}



process FASTQC {
    tag "FASTQC $srr"
    input:
    tuple val(srr), path(fastqs)
    output:
    path "$srr*zip"
    script:
    """
    $params.which_fastqc -t 5 -q $fastqs
    """
}

process MULTIQC {
    tag "MULTIQC"
    input:
    path res_dir
    // val res_dir_full
    path "*"
    output:
    val finished
    script:
    finished = true
    """
    mkdir -p \$(readlink -f ${res_dir})
    $params.which_multiqc -o ${res_dir} .
    """
}


workflow FASTQC_PIPELINE {
    take:
    fastq_paires_ch
    out_dir
    main:
    MULTIQC(out_dir,FASTQC(fastq_paires_ch).collect())
    emit:
    MULTIQC.out
}

def helpMessage() {
    log.info '''\

    
            run FASTQC-MULTiQC
            ===================================
            set fastq in dir with --indir "path/to"
            output will be writen to $params.multiqc_results_dir -
            
            further explantion:
            mandatory pramters for this pipeline is params.indir of fastq files
            params.indir=''
            recomended to set the prameters - out directory and if singleEnd:
            params.singleEnd = false
            optional prameters depend on your files structure:
            params.fastq_suffix='fastq'
            params.fastq_pat="${params.indir}/*_{1,2}.${params.fastq_suffix}"
            extra parameters:
            '''
            .stripIndent()

}

workflow {
    if (params.help){
        helpMessage()
    }
    else {  
        fastq_ch=Channel.fromFilePairs( params.fastq_pat, size: params.singleEnd ? 1 : 2 )
        .ifEmpty { System.exit(1), "Cannot find any reads matching: ${params.fastq_pat}\nIf this is single-end data, please specify --singleEnd on the command line." }
        FASTQC_PIPELINE(fastq_ch,params.multiqc_results_dir)
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
