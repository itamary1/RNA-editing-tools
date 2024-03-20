nextflow.enable.dsl=2
/*
run for example
nextflow run ./subpipelines/run_fastp.nf --indir /home/alu/twerski/Scripts/Nextflow/training_data/ggal --fastqs_outdir /home/alu/twerski/Scripts/Nextflow/training_data/two_files/fastqdedup --fastq_suffix fq
*/

process make_output_dirs {
    input:
        val prev_finished
        path out_QC_dir
        path out_fq_dir
    output:
        val finished
    script:
    finished = true
    """
    mkdir -p \$(readlink -f $out_QC_dir)
    mkdir -p \$(readlink -f $out_fq_dir)
    """

}

process FASTP {
    tag "fastp $acc"
    // for preventing overuse of memory - every fastp process using 25GB - u need 150GB 
    maxForks 6
    input:
        val create_dirs_fin
        path out_QC_dir
        path out_fq_dir
        tuple val(acc), path(fastqs)
    output:
        tuple val(acc), path('fastqs/*')
    script:
    def single = fastqs instanceof Path
    def read1 = !single ? /"${fastqs[0]}"/ : /"${fastqs}"/
    def read2 = !single ? /"${fastqs[1]}"/ : ''
    """
    r1=$read1
    r2=$read2
    mkdir fastqs
    if [ -z "\$r2" ]; then
        $params.which_fastp -i $read1 -o ${out_fq_dir}/\${r1##*/} -h ${out_QC_dir}/${acc}.html -j ${out_QC_dir}/${acc}.json $params.fastp_params
    else
        $params.which_fastp -i $read1 -I $read2 -o ${out_fq_dir}/\${r1##*/} -O ${out_fq_dir}/\${r2##*/} -h ${out_QC_dir}/${acc}.html -j ${out_QC_dir}/${acc}.json $params.fastp_params
    fi
    find -L "\$(pwd)/$out_fq_dir" -name \"${acc}*\" -type f -exec ln -s {} ./fastqs  ';'
    """
}


def helpMessage() {
    log.info '''\

    
            RUN FASTP
            ===================================
            set fastq in dir with --indir "path/to" or params.indir=
            output will be writen to params.fastp_fastqs_outdir="${launchDir}/fastq_dedup/"
            // determind where the fastp-QC file will be
            // !!!a relative!! path from outdir to fastqc
            params.outdir_to_qc="/../../Results/QUALITY_CHECK/Fastp/"
            
            >> other params: 
            params.singleEnd = false
            // prarams for finding the input fastq
            params.fastq_suffix='fastq'
            params.mates_patt='_{1,2}'
            params.fastq_pat="${params.indir}/*${params.mates_patt}.${params.fastq_suffix}"
 
            there are two config -
            1. old with dedup (no quality/length filtering)
            2. new with quality/length filtering (no dedup)

            // fastp params - >>>>>>>>>> old config
            params.fastp_dedup_params='--dedup --dup_calc_accuracy 6'
            params.fastp_disable_filtering_params='--disable_quality_filtering --disable_length_filtering --disable_trim_poly_g --disable_adapter_trimming'
            params.fastp_extra_params=''
            params.fastp_params = "$params.fastp_dedup_params $params.fastp_disable_filtering_params $params.fastp_extra_params"
           
            // fastp params - >>>>>>>>>> new config
            // if you set fastp_filter_to_length=true you must set your read data read length -
            // note that in case of fastp_filter_to_length=fasle, it will filter reads<15bp, you can disable it with --disable_length_filtering
            fastp_filter_to_length = false
            read_length = ''
            // if filter to length was set to true:
            // length trimming: trim reads with len > (params.read_length+5), to (params.read_length+5)
            fastp_length_params = params.fastp_filter_to_length ? " --max_len1 \$(($params.read_length+5)) --length_required \$(($params.read_length-10))" : ''
            // in defualt there will be no deduplication
            fastp_dedup = false
            fastp_dedup_params = params.fastp_dedup ? ' --dedup --dup_calc_accuracy 6' : ''
            // if dont dedup so dont eval to save time
            dont_eval = params.fastp_dedup ? "" : "--dont_eval_duplication"
            // > quality filtering
            filter_quality = true
            dont_filter_quality = params.filter_quality ? "" : "--dont_filter_quality"
            // quality filtering: N bases per read <= 5
            // quality filtering: average quality per read >= 30
            // quality filtering: % of low-quality bases per read <= 20
            // quality filtering: low quality per base <= 25
            // length filtering: reads len >= params.read_length-10
            fastp_quality_params= params.filter_quality ? " -n 5 -q 25 -u 20 -e 30" : ""
            fastp_extra_params= ''
            
            fastp_disable_params = " --disable_trim_poly_g --disable_adapter_trimming $params.dont_eval $dont_filter_quality"
            fastp_params = "$params.fastp_quality_params $params.fastp_dedup_params $params.fastp_disable_params $params.fastp_extra_params"
            '''
            .stripIndent()

}



workflow FASTP_PIPELINE {
    take:
    fastqCH
    QC_out_dir
    fq_outD
    prev_finished
    main:
    // we need collect for having value channel
    create_dirs_finished = make_output_dirs(prev_finished,QC_out_dir,fq_outD).collect()
    FASTP(create_dirs_finished,QC_out_dir,fq_outD, fastqCH)
    emit:
    FASTP.out
}

// workflow {
//         fastq_ch=Channel.fromFilePairs( params.fastq_pat, size: params.singleEnd ? 1 : 2 )
//         .ifEmpty { System.exit(1), "Cannot find any reads matching: ${params.fastq_pat}\nIf this is single-end data, please specify --singleEnd on the command line." }
//         FASTP_PIPELINE(fastq_ch, params.fastp_fastqs_outdir)
//         fastp_fq=FASTP_PIPELINE[0]
//         fastp_isf=FASTP_PIPELINE[1]
//         fastp_isf.view()
// }


workflow {
    if(params.help){
        helpMessage()
    } else {
    fastq_ch=Channel.fromFilePairs( params.fastq_pat, size: params.singleEnd ? 1 : 2 )
    .ifEmpty { System.exit(1), "Cannot find any reads matching: ${params.fastq_pat}\nIf this is single-end data, please specify --singleEnd on the command line." }
    FASTP_PIPELINE(fastq_ch, true).view()
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