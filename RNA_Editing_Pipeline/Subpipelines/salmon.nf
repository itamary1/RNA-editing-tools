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



def helpMessage() {
    log.info '''\
            run salmon
            ===================================
            
            set fastq dir with --indir "path/to"
            run with -profile hg38 for compatible resources, I will be glad if you add other files to othe profiles in the config
            output will be writen to params.salmon_outdir="${launchDir}/SALMON"

            params.run_salmon_summary - defualts: true - will run roni's summary script()
            (in case of docker run it will run it localy and not in docker)

            
            >>>other params:
            params.fastq_suffix='fastq'
            params.mates_patt="*_{1,2}"
            params.fastq_pat="${params.indir}/${params.mates_patt}.${params.fastq_suffix}"
            //prams of salmon index
            params.salmon_lib_type = 'A'
            params.transcripts_index = '/private/dropbox/Salmon_1.4.0/salmon_index/hg38'
            params.tx2id_geneMap = '/private/dropbox/Salmon_1.4.0/salmon_index/hg38/gencode_v32.transcriptToGeneID.tab'

            //other run params
            params.keep_other_files = false
            params.salmon_threads = 5 
            params.extra_params = ''
            '''
            .stripIndent()

}

// if(params.docker_run && params.run_salmon_summary){
//     println "you cant run salmon summary in docker! set run_salmon_summary=false"
//     System.exit(1)
// }

process CREATE_OUT_DIR {
    tag "create salmon outdir"
    input:
        path outdir
        val go
    output:
        val finished
    script:
    finished = true
    """
    mkdir -p "\$(readlink -f ${outdir})/ResultFiles/"
    mkdir -p "\$(readlink -f ${outdir})/summary/"
    """
}

process SALMON_QUANT {
    tag "salmon_quant $acc"
    maxForks params.salmon_max_workers
    input:
    val create_dir_finished
    tuple val(acc), path(fastqs)
    path salmon_transcripts_index
    path salmon_tx2id_geneMap
    path salmon_outdir
    output:
    path "salmon_result"
    script:
    def single = fastqs instanceof Path
    def read1 = !single ? /"${fastqs[0]}"/ : /"${fastqs}"/
    def read2 = !single ? /"${fastqs[1]}"/ : ''
    """
    acc_dir="\$(readlink -f ${salmon_outdir})/$acc"
    mkdir -p \$acc_dir
    if [ $single = true ]; then
        $params.which_salmon quant -l $params.salmon_lib_type -r $read1 -o \$acc_dir -i $params.transcripts_index -g $params.tx2id_geneMap -p $params.salmon_threads $params.extra_params
    else
        $params.which_salmon quant -l $params.salmon_lib_type -1 $read1 -2 $read2 -o \$acc_dir -i $params.transcripts_index -g $params.tx2id_geneMap -p $params.salmon_threads $params.extra_params
    fi
    if [ $params.move_out_salmon_files = true ]; then
        mv \$acc_dir/quant.sf \$acc_dir/../ResultFiles/${acc}.quant.sf 
        mv \$acc_dir/quant.genes.sf \$acc_dir/../ResultFiles/${acc}.quant.genes.sf
        if [ $params.keep_other_files = false ]; then
            rm -r \$acc_dir
        fi
    elif [ $params.run_salmon_summary = true ]; then
        cp \$acc_dir/quant.sf \$acc_dir/../ResultFiles/${acc}.quant.sf 
        cp \$acc_dir/quant.genes.sf \$acc_dir/../ResultFiles/${acc}.quant.genes.sf
    fi
    ln -s "\$(readlink -f ${salmon_outdir})" ./salmon_result
    """
}


process SALMON_SUMMARY {
    tag "summary salmon"
    input:
    path indir
    val prev_finished
    path salmon_res_dir
    output:
    path "${salmon_res_dir}/summary"
    script:
    """
    mkdir -p "\$(readlink -f ${salmon_res_dir})/summary" 
    $params.which_salmon_summary -i ${indir} -o ${salmon_res_dir}/summary --wide_genes --wide_samples
    # ln -s "\$(readlink -f ${salmon_res_dir})/summary" ./salmon_summary
    """
}

workflow SALMON_PIPELINE {
    take:
    fastqCH
    salmon_outdir
    ready_to_go
    main:
    salmon_results=SALMON_QUANT(CREATE_OUT_DIR(salmon_outdir,ready_to_go),fastqCH,params.transcripts_index,params.tx2id_geneMap,salmon_outdir)
    if(params.run_salmon_summary)
        SALMON_SUMMARY(salmon_results.first(),salmon_results.collect(),salmon_outdir)
    emit:
    salmon_results
}

workflow {
    if (params.help) {
        helpMessage()
    } else {
        fastq_ch=Channel.fromFilePairs( params.fastq_pat, size: params.singleEnd ? 1 : 2 ).ifEmpty { System.exit(1), "Cannot find any reads matching: ${params.fastq_pat}\nIf this is single-end data, please specify --singleEnd on the command line." }
        SALMON_PIPELINE(fastq_ch,params.salmon_outdir,true).view()
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