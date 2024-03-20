nextflow.enable.dsl=2

// >>>>>>>>>>>>>>>> star parameters  were take from /private/common/Software/PipelineConfigs/runSTAR.conf <<<<<<<<<<<<<<
params.help=false
params.have_base_config=false
if((!params.have_base_config) && (!params.help)){
 println "please attach a base config file"
 System.exit(1) 
}



File file = new File("$params.star_genome")
// println file.exists()
if(!file.exists() && !params.help) {
    //if star genome wasnt set by the user
    println "genome not exist at - $params.star_genome . make sure you choose profile and set genome_length with --genome_length, or set path genome dir with --star_genome <path/to>"
    System.exit(1)
}




def helpMessage() {
    log.info '''\
            <>>>>>>>>>>>run star<<<<<<<<<<<<<< 
            ===================================
            !!!!! when runing in a docker its asuumed that all the fastq is in the same dir !!!!
            will sort and index bam files
            will clean genome from shared memory (will wait if other process is using it)
            will run star stats to ${star_dir}/summary


            set fastq dir with --indir "path/to"
            output will be writen to params.star_outdir="${launchDir}/STAR/"

            //genome
            params.annotation='hg38' - also support mm10
            // User must select genome length!!! 
            params.genome_length='' - succsion depends of what indexes we have in our storages

            Note that clean_genome process can raise error in case of preloaded genome by other user,
             this wont affect the pipeline results

            
            >>>other params:
            // input fastq params:
            params.fastq_suffix='fastq'
            params.fastq_pat="${params.indir}/*_{1,2}.${params.fastq_suffix}"
            // if to compress unmapped
            params.compress_unmapped=true
            //if to clean star genome at the end of all proccess
            params.clean_genome=true


            // path to programs



            '''
            .stripIndent()

}


process STAR_LOAD {
    tag "load star genome"
    output:
    stdout
    script:
    """
    echo 'loading genome'
    $params.which_star --genomeDir $params.star_genome  --genomeLoad LoadAndExit
    echo 'finished loading'
    """
}
        
        
        
process STAR_alignment { 
    tag "STAR $acc"
    input:
    tuple val(acc), path(fastqs)
    val(genome_finished)
    path bams_outdir
    output:
    tuple val(acc), env(acc_outdir), path("./${acc}")
    script:
    """
    acc_outdir="\$(readlink -f ${bams_outdir})/$acc"
    mkdir -p \$acc_outdir
    $params.which_star --alignSJoverhangMin 8 --alignIntronMax 1000000 --alignMatesGapMax 600000 --outFilterMismatchNoverLmax 0.3 --outFilterMismatchNoverReadLmax 1 --outFilterMatchNminOverLread  0.66 --outFilterMultimapNmax 1 --outReadsUnmapped Fastx --outSAMattributes All --outSAMtype BAM Unsorted --quantMode GeneCounts --genomeLoad LoadAndKeep --runThreadN $params.runThreadN --genomeDir $params.star_genome --readFilesIn $fastqs --outFileNamePrefix \${acc_outdir}/${acc}.
    ln -s \$acc_outdir ./${acc}
    """
}



// for alowing only one load of the genome in the docker
// a one process for both loading and runing alignment
// in addition - all samples will be align in the same process
// TODO move to python?
process STAR_load_and_alignment_all {
    maxForks 1
    tag "load and align all"
    input:
        val ready_to_go
        // we need this input for mounting of the fastq's dir
        tuple val(acc), path(fastqs)
        val all_samples
        path bams_outdir
        // star docker wasnt able to make readlink
        val bams_outdir_full
        path genomeS
    output:
        path 'local_acc_dirs/*'
    script:
    def num_samples = all_samples.size()/2 
    // I am spliting the list to its parts - acc,read1 and read2 if paired
    // because of that i will move in the samples list in a step at the according size - 3  for paried and 2 for SE
    def step = params.singleEnd ? 2 : 3
    """
    echo 'loading genome'
    $params.which_star --genomeDir $genomeS  --genomeLoad LoadAndExit
    echo 'finished loading'
    mkdir local_acc_dirs
    ass_array=(\$(echo $all_samples | tr -d '[],'))
    echo 'star aligning $num_samples files'
    for ((i=0; i<$num_samples; i+=$params.STAR_MAX_PARALLEL)); do
        for ((j=i; j<i+$params.STAR_MAX_PARALLEL && j<$num_samples; j++)); do
            acc=\${ass_array[j*$step]}
            echo "acc \${acc}"
            ls \${ass_array[j*$step+1]}
            acc_outdir=${bams_outdir_full}/\${acc}
            mkdir -p \$acc_outdir
            ln -s \$acc_outdir local_acc_dirs/\${acc}
            # make read2 empty for case of single-end
            if [ $params.singleEnd = true ]; then
                read2=''
            else
                read2=\${ass_array[j*$step+2]}
            fi
            $params.which_star --alignSJoverhangMin ${params.alignSJoverhangMin} --alignIntronMax ${params.alignIntronMax} --alignMatesGapMax ${params.alignMatesGapMax} --outFilterMismatchNoverLmax ${params.outFilterMismatchNoverLmax} --outFilterMismatchNoverReadLmax ${params.outFilterMismatchNoverReadLmax} --outFilterMatchNminOverLread ${params.outFilterMatchNminOverLread} --outFilterMultimapNmax ${params.outFilterMultimapNmax} --outReadsUnmapped Fastx --outSAMattributes All --outSAMtype BAM Unsorted --quantMode GeneCounts --genomeLoad LoadAndKeep --runThreadN $params.runThreadN --genomeDir $genomeS --readFilesIn \${ass_array[j*$step+1]} \$read2 --outFileNamePrefix \${acc_outdir}/\${acc}. &> run_\${acc} &
        done
        wait
    done
    $params.which_star --genomeDir $genomeS  --genomeLoad Remove
    """
}

// will create channel for SORT from STAR_load_and_alignment_all output
process CH_FROM_STAR_load_and_alignment_all {
      input:
        path local_acc_dir
        path bams_outdir
      output:
        tuple val(acc), env(acc_dir), path(local_acc_dir)
      script:
        //convert local nextflow dir path object to string
        as_str="$local_acc_dir"
        // get accession  number
        acc=as_str.lastIndexOf('/').with {it != -1 ? as_str[it..<as_str.length()] : as_str}
        // get star real output dir (the local dir is just symlink to it)
        """
        acc_dir="\$(readlink -f ${bams_outdir})/$acc"
        """
}




// sort bam file by coordinates
process SORT {
    tag "SORT $acc"
    input:
    tuple val(acc), val(acc_dir), path(local_star_dir)
    path out_dir
    output:
    tuple val(acc), val(acc_dir), path(local_star_dir)
    """
    $params.which_samtools sort -l 9 -m 1G -@ 7 -o ${acc_dir}/${acc}.Aligned.sortedByCoord.out.bam ${acc_dir}/${acc}.Aligned.out.bam
    rm ${acc_dir}/${acc}.Aligned.out.bam
    """
}

// create index file for the bam (bai file)
process INDEX {
    tag "index $acc"
    input:
    tuple val(acc), val(acc_dir), path(local_star_dir)
    path out_dir
    output:
        tuple val(acc), path(local_star_dir)
    // val finished
    script:
    // finished = true
    """
    $params.which_samtools index ${acc_dir}/${acc}.Aligned.sortedByCoord.out.bam
    """
}

// compressed 'unmapped' files
process CompressUnmapped {
    tag "Compress Unmapped $acc"
    maxForks 2
    input:
    tuple val(acc), val(acc_dir), path(local_star_dir)
    path bams_outdir
    // output:
    // tuple val(acc), val(acc_dir), path(local_star_dir)
    """
    find -L $acc_dir -name '*Unmapped.out.mate?' -or -name '*ReadsPerGene.out.tab' | xargs -I magicstring bash -c 'chmod 775 magicstring'
    $params.compress_command ${acc_dir}/${acc}.ReadsPerGene.out.tab
    find -L $acc_dir -name '*Unmapped.out.mate?' | xargs -I magicstring $params.compress_command magicstring
    """
}

// remove star remain temporery files
process RemoveTemp {
    tag "RemoveTemp $acc"
    input:
    tuple val(acc), val(acc_dir), path(local_star_dir)
    path bams_outdir
    """
    rm -r ${acc_dir}/${acc}._STARtmp
    """   
}

// clean the genome from the shared memory - will wait to the last process using it
process clean_genome {
    tag "clean_genome"
    // can raise error when other use loaded the genome and we cand remove it
    errorStrategy = 'ignore' 
    input:
    val all_star_finished
    path starG
    when:
    !params.docker_run
    """
    $params.which_star --genomeLoad Remove --genomeDir $params.star_genome
    """
}

// run roni's star stats script 
process STARstats {
    tag "STARstats $bams_outdir"
    input:
        val(finished)
        path bams_outdir
    script:
        """
        $params.which_star_summarize_script -r ${bams_outdir} -o ${bams_outdir}/summary --plot -po ${bams_outdir}/summary
        """
}

workflow STAR_PIPELINE {
    take:
    ready_to_go
    fastq_ch
    bams_outdir
    main:
    // if run inside docker - run all sample alignments in one step for utilize shared memory
    if(params.docker_run){
        if(params.batch_run){
            star_dirs = STAR_alignment(fastq_ch ,STAR_LOAD().collect(),bams_outdir)
        } else {
            // align all samples and convert output to bams-channel for the next steps
            first_fastq=fastq_ch.first()
            raw_bams_dir = STAR_load_and_alignment_all(ready_to_go,first_fastq,fastq_ch.collect(),bams_outdir,bams_outdir,params.star_genome).flatten()
            star_dirs = CH_FROM_STAR_load_and_alignment_all(raw_bams_dir,bams_outdir)
        }
    } else {
        println "run star localy is not supported yet you need to add clean genome process"
        System.exit(1)
        star_dirs = STAR_alignment(fastq_ch ,STAR_LOAD().collect(),bams_outdir)
        // TODO clean genom (only in local running)
    }
    // sort and index bams
    final_bams = INDEX(SORT(star_dirs,bams_outdir),bams_outdir)
    // if flag of run stats is true
    if(params.run_star_summarize_script) 
        STARstats(final_bams.collect(),bams_outdir)
    if(params.compress_unmapped) {
        CompressUnmapped(star_dirs,bams_outdir)
        RemoveTemp(star_dirs,bams_outdir)
    }
    if(params.clean_genome) {
        clean_genome(star_dirs.collect(),params.star_genome)
    }
    emit:
    final_bams // contains tuple val(acc), path(local_star_dir)
}



workflow {
    if (params.help) {
        helpMessage()
    } else {
        fastq_ch=Channel.fromFilePairs( params.fastq_pat, size: params.singleEnd ? 1 : 2 ).ifEmpty { System.exit(1), "Cannot find any reads matching: ${params.fastq_pat}\nIf this is single-end data, please specify --singleEnd on the command line." }
        STAR_PIPELINE(true,fastq_ch,params.star_outdir)
        STAR_PIPELINE.out.final_bams.view()
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