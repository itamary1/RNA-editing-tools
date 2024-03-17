nextflow.enable.dsl=2
/*
run for example
nextflow -c /home/alu/twerski/Scripts/Nextflow/subpipelines/salmonTE_nextflow.config run  /home/alu/twerski/Scripts/Nextflow/subpipelines/salmonTE.nf --indir /home/alu/twerski/Scripts/Nextflow/training_data/real_4_files/ --salmonTE_outdir /home/alu/twerski/Scripts/Nextflow/training_data/trys/SALMPONTE/salmonTE/
*/


params.help=false

def helpMessage() {
    log.info '''\
            run salmonTE on directory with fastq - will results TPM for every transposable element
            ===================================
            this pipeline supporting run with docker only, please attach salmonTE_nextflow.config
            
            params:
            // input fastq dir
            indir=''
            //  where the results will be
            salmonTE_outdir="${launchDir}/salmonTE/"
            // quant params
            STE_threads='20'
            STE_organism='hs' // hs/mm/dm/dr ->././Danio rerio/Drosophila melanogaster
            salmonTE_extra_params=''
            salmonTE_quant_params = "--reference=$params.STE_organism --num_threads=${params.STE_threads} $params.salmonTE_extra_params" 

            // !!!!!!!!!!!!!! test step is not working on this docker, somthing broken!!!!!!!!!!!!!!!!!
            test_step = false
            conditions_file =''
            // test params:
            tabletype='csv'
            figtype='pdf'
            analysis_type='DE'
            conditions='control,treatment'
            '''
            .stripIndent()

}

params.have_base_config=false
if((!params.have_base_config) && (!params.help)){
 println "please attach a base config file"
 exit 1 
}


params.SalmonTE_config = false

if((!params.SalmonTE_config) && (!params.help)){
    println "please attach salmonTE config file using nextflow's \'-c\' flag - subpipelines/salmonTE_nextflow.config"
    exit 1
}


// !!!!!!!!!!!!!! test step is not working on this docker, somthing broken!!!!!!!!!!!!!!!!!
// salmon test == statistics
if(params.test_step) {
    println "docker test is broken"
    exit 1
}

// make sure conditions_file was set when runnig test step
if(params.test_step && !(params.conditions_file)) {
    println "please set salmonTE conditions file with  --conditions_file path/to"
    exit 1
}





// salmon test are broken in that docker
// process SALMONTE_TEST {
//     tag "salmonTE TEST(statistics test) $params.indir"
//     input:
//     val my_user
//     path('./salmonTE')
//     path "condition.csv"
//     output:
//     path('./salmonTE')
//     """
//     mkdir -p ./salmonTE/Statistics
//     mv "condition.csv" salmonTE/Statistics/
//     rm -f ./salmonTE/condition.csv
//     echo -e  packageurl=\\\"https://cran.microsoft.com/snapshot/2019-12-24/src/contrib/tidyverse_1.3.0.tar.gz\\\" > RS.R
//     echo -e "install.packages(packageurl, repos=NULL, type=\\\"source\\\")" >> RS.R
//     Rscript RS.R
//     SalmonTE.py test --inpath=SalmonTE_output --outpath=SalmonTE_statistical_test --tabletype=csv --figtype=$params.figtype --analysis_type=DE --conditions=control,treatment

//     SalmonTE.py test --inpath=salmonTE --outpath=salmonTE/Statistics --tabletype=$params.tabletype --figtype=$params.figtype --analysis_type=$params.analysis_type --conditions=$params.conditions %> salmonTE_test.out.txt
//     chown -R -L $my_user ./
//     """

// }


process SALMONTE_QUANT {
    tag "salmonTE quant $fastq_folder"
    maxForks 1
    input:
    val prev_finished
    path 'fastq_folder'
    path 'results_dir'
    output:
    path './salmonTE'
    script:
    """
    ls -l \$(readlink -f fastq_folder)
    mkdir -p \$(readlink -f results_dir)
    SalmonTE.py quant $params.salmonTE_quant_params --outpath="\$(readlink -f results_dir)" \$(readlink -f fastq_folder) &> ./salmonTE_out.txt
    ln -s \$(readlink -f results_dir) ./salmonTE
    """
}


// process MV_RESULTS {
//     tag "mv salmonTE results to $results_dir"
//     label 'simple_bash' 
//     input:
//     path 'salmonTE'
//     path 'results_dir'
//     output:
//     path('./salmonTE_local')
//     script:
//     """
//     # mkdir -p \$(readlink -f results_dir)
//     mv \$(readlink -f ./salmonTE) \$(readlink -f results_dir)
//     ln -s \$(readlink -f results_dir) ./salmonTE_local
//     """
// }


workflow SALMONTE_PIPELINE {
    take:
    prev_finished
    fastq_folderCH
    results_dir
    main:
    quant_results = SALMONTE_QUANT(prev_finished,fastq_folderCH,results_dir)
    // if(params.test_step) {
    //     condition_fileCH=Channel.fromPath(params.conditions_file)
    //     .ifEmpty { exit 1, "Cannot find your conditions_file"}
    //     MV_RESULTS(SALMONTE_TEST(the_user,quant_results,condition_fileCH))
    // } else {
    // MV_RESULTS(quant_results,results_dir)
    // }
    emit:
    SALMONTE_QUANT.out
}

workflow {
        fastq_folder=Channel.fromPath(params.indir, type:'dir')
        .ifEmpty { exit 1, "Cannot find your fastq folder"}
        SALMONTE_PIPELINE(fastq_folder,params.salmonTE_outdir).view()
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