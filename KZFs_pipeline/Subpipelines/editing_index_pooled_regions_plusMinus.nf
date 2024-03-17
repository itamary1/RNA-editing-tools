nextflow.enable.dsl=2


params.help=false


// if((!(params.one_regions_file && params.one_regions_file && params.one_RID_file)) && !params.help){
//     println "please set one regions_file and one refseq_file"
//     exit 1
// }


def helpMessage() {
    log.info '''\

            runing command example:
            
            # set vars:
            bed6_file="/home/alu/twerski/Scripts/Nextflow/Special_pipelines/Resources/KZFs/human_KZFs_Imbeault_Gencode_407Genes_CDS.bed6.bed"
            bams_dir="/private8/Projects/itamar/ZNF/jose/Differentiated_PA1/Raw_data/STAR/"
            # run:
            nohup nextflow -bg -c /home/alu/twerski/Scripts/Nextflow/Special_pipelines/Configs/Dockers/SubP_configs/editing_index_pooled_regions_plusMinus.nf.docker.config run /home/alu/twerski/Scripts/Nextflow/Special_pipelines/Subpipelines/editing_index_pooled_regions_plusMinus.nf --bams_dir $bams_dir \
            --cmp_regions_bed $cmp_regions_bed \
            --PMpool_bed6_file $bed6_file \
            -profile hg38 &> runPMP.txt &

            disjointed:
            bams=/private10/Projects/Itamar/check_Lab_pipline/try_FUP/Raw_data/STAR
            OD=/private10/Projects/Itamar/check_Lab_pipline/Try_DisJ
            script=/home/alu/twerski/Scripts/Nextflow/Special_pipelines/Subpipelines/editing_index_pooled_regions_plusMinus.nf
            conf=/home/alu/twerski/Scripts/Nextflow/Special_pipelines/Configs/Dockers/SubP_configs/editing_index_pooled_regions_plusMinus.nf.docker.config
            bed6F=/home/alu/twerski/Scripts/Nextflow/Special_pipelines/Resources/KZFs/human_KZFs_Imbeault_Gencode_407Genes_CDS_merged.bed6.sorted.bed

            nextflow -c "$conf" run "$script" --PMpool_bed6_file $bed6F --bams_dir "$bams" --PMpooledEI_result_dir "$OD" --run_disjoint -profile hg38 > run_nf.out.txt
            
            ===================================================================

            output will be writen to ${params.proj_resDir}/${regions_without_suffix} 
            - regions_without_suffix is your regions file name
            the results will be in two folders plus/minus
            you can set output dir with: both_strands_resDir
            

            plusMinusWords=["plus", "minus"]

            Eindex_result_dir="${launchDir}/Editing_index/"
            indir=''
            // !!! remeber to set your profile for othe params for editing_index.nf.docker.config
            // index run params: 
            PMkeep_cmpileup=false
            PMgenome_file=''
            PMbams_suffix='.Aligned.sortedByCoord.out.bam'
            PMper_region_output=false
            PMper_sample_output=false
            PMsnp_file=''
            PMexpression_file=''
            // index will run only on path that contain that string -
            PMmust_contain=''

            '''
            .stripIndent()

}



include { POOLED_EI_Dynamic_regions_PIPELINE as PLUS_POOLED_EI; POOLED_EI_Dynamic_regions_PIPELINE as MINUS_POOLED_EI } from "./editing_index_pooled_regions.nf"
include { POOLED_EI_Dynamic_regions_PIPELINE as POOLED_EI } from "./editing_index_pooled_regions.nf"
include { disjoint_bed_file } from "./disjoint_bed.nf"

include { COMBINE_WAITING } from "${params.levanon_lab_pipelines_dir}/Subpipelines/combine_waiting.nf" // take: waiting_one waiting_two emit: both_finished



process create_dirs_split_files {
    input:
        path result_dir
        val both_strand_dir
        path 'bed6file_regions.bed'
    output:
        env f_outdir, emit: plus_dir 
        env s_outdir, emit: minus_dir
        path './bed6file_regions.plus.bed', emit: plus_regions
        path './bed6file_regions.minus.bed', emit: minus_regions
    script:
    """
    f_outdir=${both_strand_dir}/plus
    s_outdir=${both_strand_dir}/minus
    mkdir -p \$f_outdir
    mkdir -p \$s_outdir
    $params.read_regions_command bed6file_regions.bed | awk 'BEGIN { FS = \"\\t\"}; \$6=="+" {print \$0}' | bedtools sort | bedtools merge -c 4,5,6 -o distinct > ./bed6file_regions.plus.bed
    $params.read_regions_command bed6file_regions.bed | awk 'BEGIN { FS = \"\\t\"}; \$6=="-" {print \$0}' | bedtools sort | bedtools merge -c 4,5,6 -o distinct > ./bed6file_regions.minus.bed
    """
}


// process MV_F_CREATED_BED6F {
//   input:
//     val readyToGo
//     path 'New_bed6_file'
//     path 'Outdir'
//   script:
//   """
//   mv \$(readlink -f ./New_bed6_file) ./Outdir/
//   """
// }

// process MV_S_CREATED_BED6F {
//   input:
//     val readyToGo
//     path 'New_bed6_file'
//     path 'Outdir'
//   script:
//   """
//   mv \$(readlink -f ./New_bed6_file) ./Outdir/
//   """
// }


workflow POOLED_Eindex_plusMinus_PIPELINE {
    take:
    bams_dir
    result_dir
    readyToGo
    main:
    if (params.run_disjoint){
        regions_bed_ch=disjoint_bed_file(params.PMpool_bed6_file)
        POOLED_EI(bams_dir,regions_bed_ch,params.PMpool_bed6_file,result_dir,readyToGo)
        all_finished = POOLED_EI.out
    } else {
        create_dirs_split_files(result_dir,result_dir,params.PMpool_bed6_file)
        plus_outdir=create_dirs_split_files.out.plus_dir
        minus_outdir=create_dirs_split_files.out.minus_dir
        plus_regions=create_dirs_split_files.out.plus_regions
        minus_regions=create_dirs_split_files.out.minus_regions
        F_finished = PLUS_POOLED_EI(bams_dir,plus_regions,plus_regions,plus_outdir, readyToGo)
        S_finished = MINUS_POOLED_EI(bams_dir,minus_regions,minus_regions,minus_outdir,F_finished)
        all_finished = COMBINE_WAITING(F_finished,S_finished)
    }
    emit:
    all_finished
}

// workflow POOLED_Eindex_plusMinus_PIPELINE {
//     take:
//     bams_dir
//     result_dir
//     readyToGo
//     main:
//     create_dirs_split_files(result_dir,result_dir,params.PMpool_bed6_file)
//     plus_outdir=create_dirs_split_files.out.plus_dir
//     minus_outdir=create_dirs_split_files.out.minus_dir
//     plus_regions=create_dirs_split_files.out.plus_regions
//     minus_regions=create_dirs_split_files.out.minus_regions
//     F_finished = PLUS_POOLED_EI(bams_dir,plus_regions,plus_outdir, readyToGo)
//     S_finished = MINUS_POOLED_EI(bams_dir,minus_regions,minus_outdir,F_finished)
//     all_finished = COMBINE_WAITING(F_finished,S_finished)
//     MV_F_CREATED_BED6F(F_finished,plus_regions,plus_outdir)
//     MV_S_CREATED_BED6F(S_finished,minus_regions,minus_outdir)    
//     emit:
//     all_finished
// }



workflow {
    if(params.help){
        helpMessage()
    } else {
        POOLED_Eindex_plusMinus_PIPELINE(params.bams_dir,params.PMpooledEI_result_dir, true).view()
    }
}