//Author - Itamar Twersky
nextflow.enable.dsl=2
params.help=false
params.have_base_config=false
if((!params.have_base_config) && (!params.help)){
 println "please attach a base config file"
 System.exit(1) 
}


// >>>>>>>>>>>>>>>>>>> parameters for fastp were taken from /home/alu/fulther/ConfigsHillel/downloadFixup.fastp.conf  <<<<<<<<<<<<<<<<<<<<<<<<<<<

def helpMessage() {
    log.info '''\

    
            DOWNLOAD-FIXUP PIPELINE => DOWNLOAD->FASTQC->FASTP->FASTQC
            ===================================
            nohup nextflow -bg -c /home/alu/twerski/Scripts/Nextflow/Levanon_lab_NEXTFLOW_PIPELINE/Configs/Docker/downloadAndPreprocess.nf.docker.config run /home/alu/twerski/Scripts/Nextflow/Levanon_lab_NEXTFLOW_PIPELINE/downloadAndPreprocess.nf --ACC_list /home/alu/twerski/Scripts/Nextflow/training_data/small_srr.RL75.Human.txt &> run.out.txt




            >>>>>>>>>>>>>>>>>>> paramters for fastp were taken from /home/alu/fulther/ConfigsHillel/downloadFixup.fastp.conf  <<<<<<<<<<<<<<<<<<<<<<<<<<<

            directory for your pipeline results, if wont set it wil be $PWD
            params.project_dir="$launchDir"

            >>>>flags regarding original fastq-source  
            flag for choosing to downlad the files from the sra,
            if set to false, you should set params.ACC_list file with accessions list 
            if set to true, you should set params.indir for finding fasta files 
            params.use_existing_fastq=false
            a file with sra acc list - each srr in a line >> use if use_existing_fastq was set to false <<
            params.ACC_list=''
            directory containing fastq files >> use if use_existing_fastq was set to true <<
            params.indir=''
            optional prameters depend to combine with indir:
            params.fastq_suffix='fastq'
            params.fastq_pat="${params.indir}/*_{1,2}.${params.fastq_suffix}"


            >>>>>general params
            params.raw_data_dir="${params.project_dir}/Raw_data"
            params.singleEnd=false
            params.results_dir="${params.project_dir}/Results"

            >>>>optional steps
            params.compress_original_fastq=true
            if to add step of trimming and run FASTQC again
            note that set to this flag to true will result with 3 fastq dirs and matching FASTQC - orig/dedup/dedup_trim
            params.add_trimStep=false
            // if you want to trim in addtion to dedups of fastp (in one step) - note to set num ouf bp to trim below
            params.trim_in_dedup=false

            choose how many bases to trim from fromt and tail
            params.trimstep_basesF=10
            params.trimstep_basesT=10


            '''
            .stripIndent()

}



// if (params.ACC_list=='' && !( params.use_existing_fastq &&  params.fastq_indir) && !params.help ) {
//         println "use_existing_fastq: $params.use_existing_fastq"
//         println "fastq_indir: $params.fastq_indir"
//         println "set ACC_list or set use_existing_fastq to true and choose indir(full of fastqs), run with --help for more details"
//         System.exit(1)
//     }

if ((params.compress_original_fastq && params.remove_original_fastq) && !params.help ) {
        println "set compress_original_fastq or remove_original_fastq but not both"
        System.exit(1)
    }


//take: sraACCfile emit: extract_sra.out
include { SRA_DOWNLOAD_PIPELINE } from './Subpipelines/sra_download.nf' 
//take: fastq_paires_ch,out_dir emit: MULTIQC.out
include { FASTQC_PIPELINE as  FASTQC_ORIG } from './Subpipelines/quality_check.nf'
//take: fastq_paires_ch emit: MULTIQC.out
include { FASTQC_PIPELINE as  FASTQC_PREPROCESSED } from './Subpipelines/quality_check.nf'
//take:fastqCH , prev_finished   emit: FASTP.out
include { FASTP_PIPELINE } from './Subpipelines/fastp.nf'
//take: fastq_ch,finished, emit: all_finished
include { COMPRESS_FASTQs } from './Subpipelines/compress_fastqs.nf'
//take:waiting_one waiting_two emit: result
include { COMBINE_WAITING } from './Subpipelines/combine_waiting.nf'

process make_project_dirs {
  input:
    path project_dir
  output:
    val finished
  script:
  finished = true
  """
  mkdir -p ${params.QC_dir}
  mkdir -p ${params.raw_data_dir}
  sleep 3
  """
}

process rm_fastqs {
  input:
    tuple val(acc), path(fastqs)
    val finished
  output:
    val rm_finished
  script:
    rm_finished =true
    def single = fastqs instanceof Path
    def read1 = !single ? /"${fastqs[0]}"/ : /"${fastqs}"/
    def read2 = !single ? /"${fastqs[1]}"/ : ''
    """
    r2=$read2
    rm \$(readlink -f ${read1})
    if [ ! -z "\$r2" ]; then
        rm \$(readlink -f ${read2})
    fi
    """
}

process  link_fastq{
  input:
    tuple val(acc), path(fastqs)
    path project_dir
    val dirs_finished
  output:
    tuple val(acc), path(fastqs)
  script:
    def single = fastqs instanceof Path
    def read1 = !single ? /"${fastqs[0]}"/ : /"${fastqs}"/
    def read2 = !single ? /"${fastqs[1]}"/ : ''
    """
    r2=$read2
    mkdir -p ${params.orig_fq_dir}
    ln -s ${read1} ${params.orig_fq_dir}
    if [ ! -z "\$r2" ]; then
        ln -s ${read2} ${params.orig_fq_dir}
    fi
    """
}

workflow DOWNLOAD_AND_PREPROCESS_PIPELINE {
    take:
    fastq_source
    main:
    //create project dirs
    dir_fin = make_project_dirs(params.project_dir)
    // create fastq channel from downloading or from existing dir
    if (params.use_existing_fastq){
        in_fastq_ch=Channel.fromFilePairs( fastq_source, size: params.singleEnd ? 1 : 2 ).ifEmpty { System.exit(1), "Cannot find any reads matching: ${fastq_source}\nIf this is single-end data, please specify --singleEnd on the command line." }
        orig_fastq_ch = link_fastq(in_fastq_ch,params.project_dir,dir_fin)
    } else {
        orig_fastq_ch = SRA_DOWNLOAD_PIPELINE(fastq_source,params.orig_fq_dir,dir_fin) 
    }
    // run fastqc on the original files
    qc_orig_finished = FASTQC_ORIG(orig_fastq_ch,"${params.QC_dir}/FastQC_orig_fq")

    // // run fastp to clean the original files
    preproc_fastq_ch = FASTP_PIPELINE(orig_fastq_ch,"${params.QC_dir}/FASTP_QC",params.preproc_fastq_dir,dir_fin)
    // // run fastqc on the preprocessed files
    FASTQC_PREPROCESSED(preproc_fastq_ch, "${params.QC_dir}/FastQC_preproc_fq") 
    // combine wating of two steps to one:
    origFQ_finished = COMBINE_WAITING(preproc_fastq_ch.collect(), qc_orig_finished)
    // if the user chose to compress the original files
    if (params.compress_original_fastq)
        COMPRESS_FASTQs(orig_fastq_ch , origFQ_finished)
    else {
        if(params.remove_original_fastq){
            rm_fastqs(orig_fastq_ch , origFQ_finished)
            }
        }
    
    emit:
    preproc_fastq_ch
}


workflow {
    if (params.help) {
        helpMessage()
    } else {
        if (params.use_existing_fastq) 
          DOWNLOAD_AND_PREPROCESS_PIPELINE(params.fastq_pat).view()
        else
          DOWNLOAD_AND_PREPROCESS_PIPELINE(params.ACC_list).view()
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