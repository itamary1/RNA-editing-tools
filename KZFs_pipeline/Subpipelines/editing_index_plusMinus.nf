nextflow.enable.dsl=2

/*******************

******************/


params.help=false



def helpMessage() {
    log.info '''\

          example command:
          BAMS="/private10/Projects/KZNF_itamar/KZFs_pipeline/Runs/Jose_len51/Raw_data/STAR/"
          Refseq_file=/home/alu/twerski/Scripts/Nextflow/Special_pipelines/Resources/KZFs/refseqsFormat_KrabGencodeLocations.bed
          bed6_file=/home/alu/twerski/Scripts/Nextflow/Special_pipelines/Resources/KZFs/KRABs_locations.gencode.bed
          nextflow -bg -c /home/alu/twerski/Scripts/Nextflow/Special_pipelines/Configs/Dockers/SubP_configs/editing_index_plusMinus.nf.docker.config run /home/alu/twerski/Scripts/Nextflow/Special_pipelines/Subpipelines/editing_index_plusMinus.nf -profile hg38 --bams_dir $BAMS --PMEI_bed6file $bed6_file &> run.out.txt &

            this is version of runing editing index for case of running in the both strands:
            As every region in regular editing index run will be tested on one strand according to refseq annotation.
            Here we allow you to seprate the editing index run to two different runs according to the strand you choose.

          if you use disjoint flag, it can handle overlapping regions on the same strand to get EI for every region. but! in that case you wont be able to get overall EI

          params:
          // for runinng alone (not as module)
          result_dir="${launchDir}/Editing_index/"
          bams_dir=''

          // if to run disjoint on the bed file
          // use this option when you have ovrellaping regions on the same strand
          run_disjoint = false

          read_regions_command = 'cat'
          read_refseq_command = 'cat'
          //EI resources
          PMEI_refseqfile =''
          PMEI_bed6file = ''
          // !!! remeber to set your profile for othe params for editing_index.nf.docker.config
          // index run params: 
          PMkeep_cmpileup=false
          PMgenome_file=''
          PMbams_suffix='.Aligned.sortedByCoord.out.bam'
          PMper_region_output=false
          PMper_sample_output=false
          PMsnp_file=''
          // index will run only on path that contain that string -
          PMmust_contain=''
          
            '''
            .stripIndent()

}

params.have_base_config=false
if((!params.have_base_config) && (!params.help)){
 println "please attach a base config file"
    System.exit(1)
}

File file = new File(params.PMEI_bed6file)
if((!file.exists()) && (!params.help)){
 println "file bed6 file for params.PMEI_bed6file not exist at $params.PMEI_bed6file"
    System.exit(1)
}




// EI_PIPELINE take: All_bams_dir_path,result_dir,star_finished emits:  RUN_EINDEX.out
include { EI_PIPELINE_dynamicRegRef as PLUS_EI} from "${params.levanon_lab_pipelines_dir}/Subpipelines/editing_index.nf" addParams(alu_index : false, expression_file : params.EmptyPath, keep_cmpileup : params.PMkeep_cmpileup, per_region_output : true, per_sample_output : true, run_on_multiRegionsFiles : params.run_disjoint)

include { EI_PIPELINE_dynamicRegRef as MINUS_EI} from "${params.levanon_lab_pipelines_dir}/Subpipelines/editing_index.nf" addParams(alu_index : false, expression_file : params.EmptyPath, keep_cmpileup : params.PMkeep_cmpileup, per_region_output : true, per_sample_output : true, run_on_multiRegionsFiles : params.run_disjoint)
    
include { disjoint_bed_file } from "./disjoint_bed.nf"

if(params.PMkeep_cmpileup && !(params.keep_PMfiles)){
  println "Note! keep_PMfiles flag will delete your cmpileup files"
}

process create_out_dirs {
  input:
    val go
    path result_dir
  output:
    env f_outdir, emit: plus_dir 
    env s_outdir, emit: minus_dir
  script:
  """
    f_outdir=\$(readlink -f $result_dir)/plus
    s_outdir=\$(readlink -f $result_dir)/minus
    mkdir -p \$f_outdir
    mkdir -p \$s_outdir
  """
}

process create_pseodo_refseq {
  input:
    path bed6_file
  output:
    path 'refseq_file.bed'
  script:
  """
  cat $bed6_file | awk 'BEGIN { FS=\"\\t\"} ;{OFS=\"\\t\"; print \$1,\$2,\$3,\"0\",\$4,\$6,\$2\",\",\$3\",\"}' > refseq_file.bed
  """
}

process split_refseq {
  input:
    path 'refseq_file.bed'
  output:
    path './refseq_file.plus.bed', emit: plus_refseq
    path './refseq_file.minus.bed', emit: minus_refseq
  script:
  """
  $params.read_refseq_command refseq_file.bed |  awk 'BEGIN { FS = \"\\t\"}; \$6=="+" {print \$0}' > ./refseq_file.plus.bed
  $params.read_refseq_command refseq_file.bed |  awk 'BEGIN { FS = \"\\t\"}; \$6=="-" {print \$0}' > ./refseq_file.minus.bed
  """
}

process split_regions {
  input:
    path 'bed6file_regions.bed'
  output:
    path './bedfile_regions.plus.bed', emit: plus_regions
    path './bedfile_regions.minus.bed', emit: minus_regions
  script:
  """
  # use bed3 as output - because EI demanding it -> no  -c 4,5,6 -o distinct
  $params.read_regions_command bed6file_regions.bed | awk 'BEGIN { FS = \"\\t\"}; \$6=="+" {print \$0}' | bedtools sort | bedtools merge > ./bedfile_regions.plus.bed
  $params.read_regions_command bed6file_regions.bed | awk 'BEGIN { FS = \"\\t\"}; \$6=="-" {print \$0}' | bedtools sort | bedtools merge > ./bedfile_regions.minus.bed
  """
}
  
process combine_tables {
  input:
    path 'plus_input_tables'
    path 'minus_input_tables'
    path output_dir
  output:
    path 'StrandDerivingCountsPerRegion.csv'
  script:
  """
  # delete all summary/StrandDerivingCountsPerRegion.csv - broken tables
  find -L ./plus_input_tables -name "StrandDerivingCountsPerRegion.csv" | grep summary/StrandDerivingCountsPerRegion.csv | xargs -I magS rm magS
  find -L ./minus_input_tables -name "StrandDerivingCountsPerRegion.csv" | grep summary/StrandDerivingCountsPerRegion.csv | xargs -I magS rm magS

  # ###collect and concatenate all tables
  mkdir -p \$(readlink -f ${output_dir})/summary
  # create header
  head -1 \$(find -L ./plus_input_tables -name "StrandDerivingCountsPerRegion.csv" | head -1) > ${output_dir}/summary/StrandDerivingCountsPerRegion.csv
  # add plus tables
  find -L ./plus_input_tables -name "StrandDerivingCountsPerRegion.csv" | xargs -I magS tail -n +2 magS >> ${output_dir}/summary/StrandDerivingCountsPerRegion.csv
  # add minus tables
  find -L ./minus_input_tables -name "StrandDerivingCountsPerRegion.csv" | xargs -I magS tail -n +2 magS >> ${output_dir}/summary/StrandDerivingCountsPerRegion.csv
  ln -s ${output_dir}/summary/StrandDerivingCountsPerRegion.csv StrandDerivingCountsPerRegion.csv
  # remove plus/minus dirs
  if [ $params.keep_PMfiles = false ]; then
    rm -r \$(readlink -f ./plus_input_tables) \$(readlink -f ./minus_input_tables)
  fi
  """
}

// calculate EI from reads count
// in disjoint will not calculate total EI because its not right 
process calculate_EI {
  input:
    path 'input_table'
    path 'output_dir'
  output:
    path 'output_dir'
  script:
  // in disjoint will not calculate total EI because its not right 
  calculate_totalEI = params.run_disjoint ? false : true
  run_total_EI_pythonFlag = calculate_totalEI ? 'True' : 'False'
  """
  ${params.which_calculateEI_command} ./input_table ${run_total_EI_pythonFlag} ./output_dir/summary/EditingIndex.csv
  if [ $calculate_totalEI = true ]; then
    ln -s ./output_dir/summary/EditingIndex.csv ./EditingIndex.csv
  fi
  """
}

workflow Eindex_plusMinus_PIPELINE {
    take:
    bams_dir
    results_dir
    readyToGo
    main:
    //  dirs
    create_out_dirs(readyToGo, results_dir)
    plus_outdir=create_out_dirs.out.plus_dir
    minus_outdir=create_out_dirs.out.minus_dir
    // split refseq
    if(params.PMEI_refseqfile) {
      refseq_file = params.PMEI_refseqfile
    } else {
      refseq_file = create_pseodo_refseq(params.PMEI_bed6file)
    }
    split_refseq(refseq_file)
    plus_refseq=split_refseq.out.plus_refseq
    minus_refseq=split_refseq.out.minus_refseq
    if (params.run_disjoint){
        disjoint_bed_file(params.PMEI_bed6file)
        disjoint_bed_file.out.flatten().map { it -> "$it"} .branch {
            plus: it.contains('plus')
            minus: it.contains('minus')
        }
        .set { bed_disjointedPM }
      plus_regions = bed_disjointedPM.plus
      minus_regions = bed_disjointedPM.minus
    } else { 
      split_regions(params.PMEI_bed6file)
      plus_regions = split_regions.out.plus_regions
      minus_regions = split_regions.out.minus_regions
    }
    PLUS_EI_res = PLUS_EI(bams_dir,plus_refseq,plus_regions,plus_outdir, readyToGo)
    MINUS_EI_res = MINUS_EI(bams_dir,minus_refseq,minus_regions,minus_outdir, PLUS_EI_res)
    final_table = combine_tables(PLUS_EI_res,MINUS_EI_res,results_dir)
    all_finished=calculate_EI(final_table,results_dir)
    emit:
    all_finished
}



workflow {
    if (params.help)
        helpMessage()
    else {
        // println "$params.bams_dir , $params.result_dir"
        Eindex_plusMinus_PIPELINE(params.bams_dir,params.result_dir, true)
        Eindex_plusMinus_PIPELINE.out.view()
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
