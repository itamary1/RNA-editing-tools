//Author - Itamar Twersky
nextflow.enable.dsl=2
params.help=false
params.have_base_config=false
if((!params.have_base_config) && (!params.help)){
 println "please attach a base config file"
 System.exit(1) 
}





def helpMessage() {
    log.info '''\
            SRR-KZFs analize pipeline => DOWNLOAD->FASTQC->FASTP->FASTQC(TRIM->FASTQC)
            ===================================                   
            >>>>>>>>>>>>>>>>>>> paramters for fastp were taken from samboseq projet TODO   <<<<<<<<<<<<<<<<<<<<<<<<<<<

             command for example:
             conf=/home/alu/twerski/Scripts/Nextflow/Special_pipelines/Configs/Dockers/analize_editing_KZF_fullPipeline.nf.docker.config
             nfScript=/home/alu/twerski/Scripts/Nextflow/Special_pipelines/analize_editing_KZF_fullPipeline.nf
             nextflow -bg -c $conf run $nfScript -profile hg38 --use_existing_fastq --fastq_indir /home/alu/twerski/Scripts/Nextflow/training_data/jose_oneSamp --genome_length 50 --project_dir $PWD  --compress_original_fastq false &> run_log.txt

            directory for your pipeline results, if wont set it wil be $PWD
            params.project_dir="$launchDir"

            >>>>flags regarding original fastq-source  
            flag for choosing to downlad the files from the sra,
            if set to true, you should set params.ACC_list file with accessions list 
            if set to false, you should set params.fastq_indir for finding fasta files 
            params.download_srrs=true
            a file with sra acc list - each srr in a line >> use if download_srrs was set to true <<
            params.ACC_list=''
            directory containing fastq files >> use if use_existing_fastq was set to true <<
            params.fastq_indir=''
            optional prameters depend to combine with fastq_indir:
            params.fastq_suffix='fastq'
            params.fastq_pat="${params.fastq_indir}/*_{1,2}.${params.fastq_suffix}"

            '''
            .stripIndent()

}

if (params.help) {
        helpMessage()
        System.exit(0)
    }

params.profile_selected=''

if((!params.profile_selected) && (!params.help)){
 println "please select profile"
 System.exit(1) 
}

include { DOWNLOAD_AND_PREPROCESS_PIPELINE } from "${params.levanon_lab_pipelines_dir}/downloadAndPreprocess.nf" // takr dir or list, emit fastqCH
include { ANALIZE_EDITING_PIPELINE } from "${params.levanon_lab_pipelines_dir}/analizeEditing_on_cleanData.nf" addParams(run_known : false, compress_fastqs : false, run_salmon_summary: true, move_out_salmon_files: false)
include {SALMONTE_PIPELINE} from "./Subpipelines/salmonTE.nf"


// >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> look on all the kzfs transcripts

// run EI on all tx for total EI
include { Eindex_plusMinus_PIPELINE as EI_total_KZFs_tx } from "./Subpipelines/editing_index_plusMinus.nf" addParams(PMEI_bed6file : params.KZF_tx_EI_bed6Regions)

// run EI disjont+multi_regions on tx for get per gene tx' EI
include { Eindex_plusMinus_PIPELINE as EI_perGene_KZFs_tx } from "./Subpipelines/editing_index_plusMinus.nf" addParams(run_disjoint : true, PMEI_bed6file : params.KZF_tx_EI_bed6Regions)

// for running editing index on all the KZFs that contains editing sites orshai list
include { Eindex_plusMinus_PIPELINE as EI_tx_OrshList_KZFs } from "./Subpipelines/editing_index_plusMinus.nf" addParams(PMEI_bed6file : params.KZF_contains_OrsahiSites_tx_regions)

// for running editing index on all the KZFs that in orshai list have more then 10 sites
include { Eindex_plusMinus_PIPELINE as EI_OrshListEnriched_KZFs_tx } from "./Subpipelines/editing_index_plusMinus.nf" addParams(PMEI_bed6file : params.KZF_enriched_OrsahiSites_tx_regions)

// for running editing index on all KZF that have inverted neighbor KZF
include { Eindex_plusMinus_PIPELINE as EI_inverted_KZFs_tx } from "./Subpipelines/editing_index_plusMinus.nf" addParams(PMEI_bed6file : params.KZF_inveterd_tx_regions)

// >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> look on the exons
// run EI on all EXONS for total EI
include { Eindex_plusMinus_PIPELINE as EI_total_KZFs_exons } from "./Subpipelines/editing_index_plusMinus.nf" addParams(PMEI_bed6file : params.KZF_exons_EI_bed6Regions)

// per_gene
include { POOLED_Eindex_plusMinus_PIPELINE as KZF_exons_pooledEI } from "./Subpipelines/editing_index_pooled_regions_plusMinus.nf" addParams(PMpool_bed6_file : params.KZF_exons_EI_bed6Regions, run_disjoint : true)

// KZFs that contains editing sites orshai list
include { Eindex_plusMinus_PIPELINE as EI_exons_OrshList_KZFs } from "./Subpipelines/editing_index_plusMinus.nf" addParams(PMEI_bed6file : params.KZF_contains_OrsahiSites_exons_regions)


// have more then 10 sites
include { Eindex_plusMinus_PIPELINE as EI_OrshListEnriched_KZFs_exons } from "./Subpipelines/editing_index_plusMinus.nf" addParams(PMEI_bed6file : params.KZF_enriched_OrsahiSites_exons_regions)

// for running editing index on all KZF that have inverted neighbor KZF
include { Eindex_plusMinus_PIPELINE as EI_inverted_KZFs_exons } from "./Subpipelines/editing_index_plusMinus.nf" addParams(PMEI_bed6file : params.KZF_inveterd_exons_regions)

// >>>>>>>>>>>>>>>>>>>>>>>> look on the CDS 

// run EI on total KZF_CDSs for total EI
include { Eindex_plusMinus_PIPELINE as EI_total_KZFs_CDS } from "./Subpipelines/editing_index_plusMinus.nf" addParams(PMEI_bed6file : params.KZF_CDS_EI_bed6Regions)

// run disjoint_pooled_EI for get per gene CDS EI
include { POOLED_Eindex_plusMinus_PIPELINE as KZF_CDS_pooledEI } from "./Subpipelines/editing_index_pooled_regions_plusMinus.nf" addParams(PMpool_bed6_file : params.KZF_CDS_EI_bed6Regions, run_disjoint : true)

// run editing index on all KZFs cds that in orshai list
include { Eindex_plusMinus_PIPELINE as EI_CDS_OrshList_KZFs } from "./Subpipelines/editing_index_plusMinus.nf" addParams(PMEI_bed6file : params.KZF_contains_OrsahiSites_CDS_regions)

// have more then 10 sites
include { Eindex_plusMinus_PIPELINE as EI_OrshListEnriched_KZFs_CDS } from "./Subpipelines/editing_index_plusMinus.nf" addParams(PMEI_bed6file : params.KZF_enriched_OrsahiSites_CDS_regions)

// for running editing index on all KZF that have inverted neighbor KZF
include { Eindex_plusMinus_PIPELINE as EI_inverted_KZFs_CDS } from "./Subpipelines/editing_index_plusMinus.nf" addParams(PMEI_bed6file : params.KZF_inveterd_CDS_regions)

// -------------------------------- Orshai sites in KZFs' CDS EditingIndex

// run EI on total KZF_CDSs for total EI
include { Eindex_plusMinus_PIPELINE as EI_total_KZFs_OrshaiInCDS } from "./Subpipelines/editing_index_plusMinus.nf" addParams(PMEI_bed6file : params.KZF_OrshaiInCDS_EI_bed6Regions)

// run disjoint_pooled_EI for get per gene CDS EI
include { POOLED_Eindex_plusMinus_PIPELINE as KZF_OrshaiInCDS_pooledEI } from "./Subpipelines/editing_index_pooled_regions_plusMinus.nf" addParams(PMpool_bed6_file : params.KZF_OrshaiInCDS_EI_bed6Regions, run_disjoint : true)

// // ------------------------------------------------> check CDS signal by random control genes
include { Eindex_plusMinus_PIPELINE as EI_CDS_GROUP_A } from "./Subpipelines/editing_index_plusMinus.nf" addParams(PMEI_bed6file : params.controlA_CDS_regions)
include { Eindex_plusMinus_PIPELINE as EI_CDS_GROUP_B } from "./Subpipelines/editing_index_plusMinus.nf" addParams(PMEI_bed6file : params.controlB_CDS_regions)
include { Eindex_plusMinus_PIPELINE as EI_CDS_GROUP_C } from "./Subpipelines/editing_index_plusMinus.nf" addParams(PMEI_bed6file : params.controlC_CDS_regions)
// try again with original EI

include { EI_PIPELINE as EI_CDS_GROUP_A_Orig } from "${params.levanon_lab_pipelines_dir}/Subpipelines/editing_index.nf" addParams(alu_index : false, keep_cmpileup : false, per_region_output : false, per_sample_output : false, regions : params.controlA_CDS_regions_bed3)







include { COMBINE_WAITING as COMB1; COMBINE_WAITING as COMB2; COMBINE_WAITING as COMB3; COMBINE_WAITING as COMB4;COMBINE_WAITING as COMB5 } from "${params.levanon_lab_pipelines_dir}/Subpipelines/combine_waiting.nf" //    take: waiting_one waiting_two emit: both_finished

include { COMPRESS_FASTQs } from "${params.levanon_lab_pipelines_dir}/Subpipelines/compress_fastqs.nf"


workflow KZF_PIPELINE {
    take:
    fastq_source
    main:
    
    if(params.debug_kzf_pipeline){
        labPipline_finished = true
        all_finished = true
    // run Levanon-lab's analize editing pipeline
    } else  {    
        clean_fastqCH=DOWNLOAD_AND_PREPROCESS_PIPELINE(fastq_source)
        ANALIZE_EDITING_PIPELINE(clean_fastqCH)
        labPipline_finished = ANALIZE_EDITING_PIPELINE.out.everything_finished
        bamsDirs_ch = ANALIZE_EDITING_PIPELINE.out.bams_dirs_ch
        fastp_finished = clean_fastqCH.collect()

    }


    bams_dir = params.star_dir
    
    // // // ---------------------------------------------------- TE analysis
    // // run salmonTE on the clean fastq dir (path from the downloadAndPreprocess config)
    // // note that salmonTE is built to work on dir of fastqs
    // SalmonTE_finished = SALMONTE_PIPELINE(fastp_finished,params.preproc_fastq_dir,params.SalmonTE_result_dir)
    // fastq_finished = COMB1(labPipline_finished,SalmonTE_finished)

    // // ----------------------------------------------------------------------------------------------------------> analize editing in KZFs

    // // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> look on all the kzfs transcripts
    // // run on all KZFs transcripts
    // // run twice:  PMsplited for looking per_gene and altogether
    // EI_txPM_finished = EI_perGene_KZFs_tx(bams_dir,params.KZF_tx_EI_perGene_dir,labPipline_finished)
    // totalEI_tx_fin = EI_total_KZFs_tx(bams_dir,params.KZF_tx_EI_dir,labPipline_finished)
    // // run on all KZFs containing sites in orshai list transcripts
    // EI_tx_OrshList_fin = EI_tx_OrshList_KZFs(bams_dir,params.KZF_contains_OrsahiSites_tx_EI_dir,EI_txPM_finished)
    // // run on all KZFs containing more then 10 sites in orshai list transcripts
    // EI_enriched_OrsahiSites_tx_fin = EI_OrshListEnriched_KZFs_tx(bams_dir,params.KZF_enriched_OrsahiSites_tx_EI_dir,totalEI_tx_fin)
    // tx_finished=COMB2(EI_enriched_OrsahiSites_tx_fin,EI_tx_OrshList_fin)

    // // run editing index on all KZF that have inverted neighbor KZF
    // EI_inverted_KZFs_tx_fin = EI_inverted_KZFs_tx(bams_dir,params.KZF_inveterd_tx_EI_dir,totalEI_tx_fin)

    // // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> look on all the kzfs exones

    // // run on all KZFs exons
    // // run twice: altogether and Pooled for looking per_gene
    // totalEI_exons_finished = EI_total_KZFs_exons(bams_dir,params.KZF_exons_EI_dir,tx_finished)
    // KZF_exons_pooledEI_finished =KZF_exons_pooledEI(bams_dir,params.KZF_exons_pooledEI_dir,tx_finished)
    // // run on all KZFs containing sites in orshai list
    // EI_exons_OrshList_KZFs_finished = EI_exons_OrshList_KZFs(bams_dir, params.KZF_contains_OrsahiSites_exons_EI_dir, totalEI_exons_finished)
    // // run on all KZFs containing more then 10 sites in orshai list
    // EI_OrshListEnriched_KZFs_exons_finished= EI_OrshListEnriched_KZFs_exons(bams_dir,params.KZF_enriched_OrsahiSites_exons_EI_dir, KZF_exons_pooledEI_finished)

    // // run on all KZF that have inverted neighbor KZF
    // EI_inverted_KZFs_exons_finished= EI_inverted_KZFs_exons(bams_dir,params.KZF_inveterd_exons_EI_dir, KZF_exons_pooledEI_finished)

    // exones_finished=COMB3(EI_exons_OrshList_KZFs_finished,EI_OrshListEnriched_KZFs_exons_finished)

    // // // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> look on all the kzfs cds
    // // // run on all KZFs CDS
    // // // run twice: altogether and Pooled for looking per_gene
    // totalEI_CDS_finished = EI_total_KZFs_CDS(bams_dir,params.KZF_CDS_EI_dir,exones_finished)
    // KZF_CDS_pooledEI_finished =KZF_CDS_pooledEI(bams_dir,params.KZF_CDS_pooledEI_dir, exones_finished)
    // // run on all KZFs containing sites in orshai 
    // EI_CDS_OrshList_KZFs_finished = EI_CDS_OrshList_KZFs(bams_dir, params.KZF_contains_OrsahiSites_CDS_EI_dir, totalEI_CDS_finished)
    // // run on all KZFs containing more then 10 sites in orshai list
    // EI_OrshListEnriched_KZFs_CDS_finished = EI_OrshListEnriched_KZFs_CDS(bams_dir, params.KZF_enriched_OrsahiSites_CDS_EI_dir, KZF_CDS_pooledEI_finished)

    // // run on all KZF that have inverted neighbor KZF
    // EI_inverted_KZFs_CDS_finished = EI_inverted_KZFs_CDS(bams_dir, params.KZF_inveterd_CDS_EI_dir, KZF_CDS_pooledEI_finished)
    

    // CDS_finihed = COMB4(EI_OrshListEnriched_KZFs_CDS_finished,EI_CDS_OrshList_KZFs_finished)
    // // Orshai sites in KZFs' CDS EditingIndex    
    // // total  EI
    // EI_total_KZFs_OrshaiInCDS_finished =EI_total_KZFs_OrshaiInCDS(bams_dir,params.KZF_OrshaiInCDS_EI_dir,CDS_finihed)
    
    // //    per_gene
    // // !!!!!!!!!!some bug in pooling scripts make this crashing
    // // KZF_OrshaiInCDS_pooledEI_finished =KZF_OrshaiInCDS_pooledEI(bams_dir,params.KZF_OrshaiInCDS_pooledEI_dir,CDS_finihed)
    // // OrshaiInCDS_finished = COMB5(EI_total_KZFs_OrshaiInCDS_finished,KZF_OrshaiInCDS_pooledEI_finished)


    // run EI on CDS of 3 random groups of genes as control
    EI_CDS_GROUP_A_f=EI_CDS_GROUP_A(bams_dir,params.controlA_CDS_EI_dir,true)
    EI_CDS_GROUP_B_f=EI_CDS_GROUP_B(bams_dir,params.controlB_CDS_EI_dir,true)
    EI_CDS_GROUP_C_f=EI_CDS_GROUP_C(bams_dir,params.controlC_CDS_EI_dir,true)
    
    // if(params.compress_fastqs_after_analysis && ! params.debug_kzf_pipeline)
    //     COMPRESS_FASTQs(clean_fastqCH, fastq_finished)
    emit:
    EI_CDS_GROUP_C_f
}

workflow {
    if (params.use_existing_fastq){
        KZF_PIPELINE(params.fastq_pat).view()
    } else {
        KZF_PIPELINE(params.ACC_list).view()
    }
}