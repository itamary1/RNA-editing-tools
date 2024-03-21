#    ----------------------------------------------------------------------------------------------------------
#
#                                   RNA Editing Standard Pipeline    
#
#   -----------------------------------------------------------------------------------------------------------

An pipeline designed to run RNA editing analysis.

the pipeline input is sra accessions list or folder with fastq files

every .nf file have compatible .config file in the Configs folder

### run the nf file with:

nextflow -c file.config run file.nf 

set your deserved genome and dependencies by :
-profile
for example: -profile hg38

if you run parts that including STAR subpipeline you also should run 
----genome_length ${length}

we supporting length of 50,75,,100,125,150


## for running the whole pipeline run(example command):
    conf=${P_DIR}/Configs/Dockers/analize_editing_KZF_fullPipeline.nf.docker.config
    nfScript=${P_DIR}/analize_editing_KZF_fullPipeline.nf
    nextflow -bg -c $conf run $nfScript -profile hg38 --use_existing_fastq --fastq_indir ${P_DIR}/Scripts/Nextflow/training_data/jose_oneSamp --genome_length 50 --project_dir $PWD  --compress_original_fastq false &> run_log.txt



## the pipeline consists of two major pipelines:
-  downloadAndPreprocess.nf
    will run quality checks and will clean the retrieved fastq files
-  analizeEditing_on_cleanData.nf
    will map fastq files using star, and calculate gene expression and RNA editing in several ways
you can run any of the two pipeline by itself or run the small subpipelines in the Subpipline folder

every nf file contain example of running command if you run its with --help


# installation:

    install Nextflow and all its dependency https://nextflow.io/ 

    run init.sh file - it will download all the heavy resources from the web and will create salmon and star indexes
    !!!!  run init.sh wil take many hours !!!!


##    If you not run on levanon labs servers:

        You can cant run the summarizing scripts, so keep the params:
        run_cmp_summarize = false
        run_salmon_summary = false
        run_star_summarize_script = false

 
