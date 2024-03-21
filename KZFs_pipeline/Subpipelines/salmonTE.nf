nextflow.enable.dsl=2


params.help=false

def helpMessage() {
    log.info '''\
            run salmonTE on directory with fastq - will results TPM for every transposable element
            ===================================
            this pipeline supporting run with docker only, please attach salmonTE_nextflow.config
            
            command example:
            
            conf=${P_DIR}/Configs/Dockers/SubP_configs/salmonTE.nf.docker.config
            scr=${P_DIR}/Subpipelines/salmonTE.nf
            indir=${P_DIR}/Scripts/Nextflow/training_data/trys/STE_test/try_buildEpim/SalmonTE/example
            outdir=${P_DIR}/temp/try_CMP_nf/Results
            nextflow -c $conf run $scr --salmonTE_outdir $outdir --indir $indir --test_step --conditions_file ${P_DIR}/temp/try_CMP_nf/condition.csv


            all params and defualts:
            // input fastq dir
            indir=''
            //  where the results will be
            salmonTE_outdir="${launchDir}/salmonTE/"
            // quant params
            STE_threads='20'
            STE_organism='hs' // hs/mm/dm/dr ->././Danio rerio/Drosophila melanogaster
            salmonTE_extra_params=''
            salmonTE_quant_params = "--reference=$params.STE_organism --num_threads=${params.STE_threads} $params.salmonTE_extra_params" 

            // if you want to run quant and test:
            test_step = false
            // to run test without quant - in that case you should give quant_dir
            test_step_only = false
            quant_dir = ''
            // test params:
            // in any case of runing test you need conditions file
            // this is CSV contaning the fields filename,condition(no header)
            // filename - without suffx, in paried end - one row for both files without file number suffix
            // you must have exactly 2 conditions
            // your first condition after sort will account as control in the DE analysis if you will not set a one
            control_condition=''
            conditions_file ='' 
            STE_tableType='csv'

            // you can set in additional
            //analysis_type=LM dosnt work duo to docker problems
            // --figtype= --analysis_type=
            STE_test_additional = ''

                    '''
            .stripIndent()

}

params.have_base_config=false
if((!params.have_base_config) && (!params.help)){
    println "please attach a base config file"
    System.exit(1)
}


params.SalmonTE_config = false

if((!params.SalmonTE_config) && (!params.help)){
    println "please attach salmonTE config file using nextflow's \'-c\' flag - subpipelines/salmonTE_nextflow.config"
    System.exit(1)
}



// make sure conditions_file was set when runnig test step
if((params.test_step || params.test_step_only) && !(params.conditions_file)) {
    println "please set salmonTE conditions file with  --conditions_file path/to"
    System.exit(1)
}


if(params.test_step && params.test_step_only) {
    println "if you want to run --test_step_only, set test_step to false and test_step_only to true."
    println "test_step is for run quant and then test. exiting.."
    System.exit(1)
}

if((params.test_step_only && !(params.quant_dir)) || (params.quant_dir && !(params.test_step_only))) {
    println "only when using --test_step_only you should(must) set --quant_dir, exiting..."
    System.exit(1)
}


process GET_USER_and_MKDIRS{
    tag "get user"
    input:
        path 'results_dir'
        val go
    output:
        env my_user
    script:
    """
    mkdir -p \$(readlink -f results_dir)
    my_user=\$(id -u \${USER}):\$(id -g \${USER})
    """
}

process SALMONTE_QUANT {
    tag "salmonTE quant $fastq_folder"
    maxForks 1
    input:
    val prev_finished
    val my_user
    path 'fastq_folder'
    path 'salmonTE_quant_result'
    output:
    path './salmonTE_quant_result'
    script:
    """
    # ls -l \$(readlink -f fastq_folder
    full_res_dir=\$(readlink -f salmonTE_quant_result)
    mkdir -p \${full_res_dir}
    /opt/SalmonTE/SalmonTE.py quant $params.salmonTE_quant_params --outpath="\${full_res_dir}" \$(readlink -f fastq_folder) &> STE.out.txt || echo STE_finished_with_error   
    echo chowning
    chown -R $my_user ./STE.out.txt && chmod -R u+rw ./STE.out.txt
    rm -r ./.snake*
    chown -R $my_user \${full_res_dir}
    chmod -R u+rw \${full_res_dir}
    """
}


process SALMONTE_TEST {
    tag "salmonTE TEST(statistics test) $params.indir"
    input:
    path 'quant_result'
    path 'salmonTE_test_results'
    path "conditions_table.csv"
    output:
    path './salmonTE_test_results'
    """
    full_res_dir=\$(readlink -f salmonTE_test_results)
    # if there is no condition file in the quant_result dir its probably not the right dir - exit
    if [ ! -f ./quant_result/condition.csv ]; then
        echo "there are not condition file in the quant_result dir, exiting.."
        exit 1
    fi
    # save orginal file
    mv ./quant_result/condition.csv ./quant_result/template_condition_file.csv

    # create updated condition file
    # select control condition
    if [ -z $params.control_condition ]; then
        control_condition=\$(tail -n +2 conditions_table.csv | cut -d "," -f2 | sort | uniq | head -1)
    else
        control_condition=$params.control_condition
    fi
    #save user file to dict
    awk 'BEGIN{FS=",";OFS=","
    while(( getline line<"conditions_table.csv") > 0 ) {
        split(line, tempArr, ",")
        condition_dict[tempArr[1]] = tempArr[2]
    }}
    # map condition to file
    NR>1 {print \$1,condition_dict[\$1]}' ./quant_result/template_condition_file.csv | awk -F "," -v control_condition="\$control_condition" 'BEGIN{OFS=","; print "SampleID,condition"};{if(\$2==control_condition) { print \$1,"control" } else { print \$1,"treatment" }}' > ./quant_result/condition.csv

    

    # get the uniq set of conditions and the number of them
    num_conditions=\$(tail -n +2 ./quant_result/condition.csv | cut -d "," -f2 | sort | uniq | wc -l)
    if ! [ \$num_conditions = 2 ]; then
        echo 'ERROR number of conditons must be 2'
        echo "there is \$num_conditions conditions"
        exit 1
    fi

    # commented because I us control,treatment conditions
    # all_conditions=\$(tail -n +2 ./quant_result/condition.csv | cut -d "," -f2 | sort | uniq  | grep -v "control")
    # all_conditions="control,\${all_conditions}" 
    all_conditions="control,treatment"

    # put the information about which is test and which is treatment in the 
    cp ./quant_result/condition.csv \${full_res_dir}/condition.csv
    echo "STE_command: /opt/SalmonTE/SalmonTE.py test --inpath=\$(readlink -f quant_result) --outpath=\$full_res_dir --conditions=\$all_conditions --tabletype=${params.STE_tableType} ${params.STE_test_additional}" > \$full_res_dir/STE_command.sh
    #
    mkdir -p \${full_res_dir}
    # run test
    /opt/SalmonTE/SalmonTE.py test --inpath=\$(readlink -f quant_result) --outpath=\$full_res_dir --conditions=\$all_conditions --tabletype=${params.STE_tableType} ${params.STE_test_additional}
    """

}


workflow SALMONTE_PIPELINE {
    take:
    prev_finished
    fastq_folderCH
    results_dir
    main:
    if (params.test_step_only){
            condition_fileCH=Channel.fromPath(params.conditions_file)
            .ifEmpty { exit 1, "Cannot find your conditions_file"}
            SALMONTE_TEST(params.quant_dir,results_dir, condition_fileCH)
            res = SALMONTE_TEST.out
    } else {
        the_user=GET_USER_and_MKDIRS(results_dir,prev_finished)
        quant_results = SALMONTE_QUANT(prev_finished,the_user,fastq_folderCH,"${results_dir}/Quant")
        if(params.test_step) {
            condition_fileCH=Channel.fromPath(params.conditions_file)
            .ifEmpty { exit 1, "Cannot find your conditions_file"}
            SALMONTE_TEST(quant_results,"${results_dir}/Test", condition_fileCH)
            res = SALMONTE_TEST.out
        } else { res=SALMONTE_QUANT.out }
    }
    emit:
    res
}

// workflow SALMONTE_GROUPS_PIPELINE {
//     split_to_groups(params.groups_metadata_csv)
//     fasta_group_ch = split_to_groups.out.fasta
//     condition_group_ch = split_to_groups.out.conditions
// }




workflow {
    if(params.test_step_only) {
        SALMONTE_PIPELINE(true,true,params.salmonTE_outdir).view()
    } else {
        fastq_folder=Channel.fromPath(params.indir, type:'dir')
        .ifEmpty { exit 1, "Cannot find your fastq folder at $params.indir"}
        SALMONTE_PIPELINE(true,fastq_folder,params.salmonTE_outdir).view()
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