nextflow.enable.dsl=2
// >>>>>>>>>>>>>>>>>>> paramters for fastp were taken from /home/alu/fulther/ConfigsHillel/downloadFixup.fastp.conf  <<<<<<<<<<<<<<<<<<<<<<<<<<<

// to run this wirkflow you only need to run with --ACC_list "path/to/file.txt" --outdir "path/to/dir" 
params.help=false
params.have_base_config=false
if((!params.have_base_config) && (!params.help)){
 println "please attach a base config file"
 System.exit(1) 
}

if(params.NGC){
    NGC_flag="--ngc $NGC"
} else {
    NGC_flag=""
}



def helpMessage() {
    log.info '''\
            

            download from sra
            ===================================
            the only params
            params.ACC_list=''
            params.outdir="$launchDir/fastq"
            !!! dont set outdir to some path uptree to $workDIr !!! - there is bug in fasterq-dump.2.10.8 
            
            when downloading paired samples, somtimes there is third file for unpaired reads - flag to delete it
            params.rm_unpairedFile=true

            extra prefetch params:
            prefetch_extra_p = ""

            // Option to give prefetch ngc file  (.ngc) for downloading from dbGaP
            // NOTE!!! I never tried this option - maybe U need to mount ngc path in case of using docker
            params.NGC=""
            '''
            .stripIndent()

}



process download_sequnetial {
    // errorStrategy 'retry'
    // maxRetries 4
    tag "Download $acc"
    maxForks 1
    input:
    val(acc)
    path 'fq_outdir'
    // val fq_outdir_full
    val(prev_finished)
    output:
    tuple val(acc), path("temp"), path('fq_outdir')
    script:
    """
    mkdir temp
    # creating the outdir here because sequential download preventing race condition
    mkdir -p \$(readlink -f fq_outdir)
    # download compressed sra file - somtimes its put sra files not inside directory so I am using roni's command
    # $params.which_prefetch -C yes $acc  -O ./temp/ $params.prefetch_extra_p $NGC_flag &> ./sra_toolkit.out.txt

    # roni command  (see sra docs for good practice of using directory name as output)
    $params.which_prefetch -C yes $acc -o ./temp/${acc}.sra $params.prefetch_extra_p $NGC_flag &> ./temp/sra_toolkit.out.txt

    # check download success
    if cat ./temp/sra_toolkit.out.txt | grep -q "successfully" ; then
        echo "download successfully"
        if cat ./temp/sra_toolkit.out.txt | grep -q "sralite" ; then
            echo -e "lite version was downloaded \n exiting.."
            exit 1
        fi
    else
        echo -e "download unsuccessfuly \n exiting.."
        exit 1
    fi
    """
}

process extract_sra {
    tag "extract $acc"

    input:
    tuple val(acc), path(sra_dir), path('fq_outdir')
    output:
    tuple val(acc), path('fastqs/*')
    """
    # extract_sra file
    echo "start extracting"
    {
        $params.which_fasterq_dump -e 6 -O fq_outdir $NGC_flag ./temp/${acc}.sra && echo "finished extracting, now deleting redundant files" && rm -r \$(readlink -f ./temp)
    } || {
        echo "somthing went wrong while extracting, exiting..."
        exit 1
    }    
    # delete unpaired file if exist
    if [ $params.rm_unpairedFile = true ] && compgen -G "fq_outdir/${acc}*_1*" > /dev/null && compgen -G "fq_outdir/${acc}*_2*" > /dev/null; then
        if [ -f fq_outdir/${acc}.fastq ]; then
            echo "removing unpaired reads file because of flag was set to true, file list before:"
            ls fq_outdir
            rm \$(readlink -f fq_outdir/${acc}.fastq) 
        fi
    fi
    # create local direcotory with fastqs links for passing the files to the next steps of the pipeline
    mkdir fastqs
    find -L "\$(pwd)/fq_outdir" -name \"${acc}*\" -type f -exec ln -s {} ./fastqs  ';'
    """
} 

/*
this pipeline take list of ["srr_acc_file", "output_dir"]
*/
workflow SRA_DOWNLOAD_PIPELINE {
    take:
    sraACCfile
    fastq_o
    previous_finished
    main:
    if (params.ACC_list) {
    srrs = Channel.fromPath(sraACCfile).splitText().map { it.replaceFirst(/\n/,'') }.ifEmpty{ System.exit(1), "unable to find accessions in path $params.ACC_list" }
    extract_sra(download_sequnetial(srrs,fastq_o,previous_finished)) 
    }
    else {
        println "\nplease set path to sra accessions file\n"
        System.exit(1)
    }
    emit:
    extract_sra.out
}


workflow {
    if (params.help) {
        helpMessage()
    } else {
        if (params.ACC_list) {
            SRA_DOWNLOAD_PIPELINE(params.ACC_list,params.sra_fastq_outdir,true).view()
            }
        else {
            println "\nplease set path to sra accessions file\n"
            System.exit(1)
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