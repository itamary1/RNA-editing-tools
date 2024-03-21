
#    ----------------------------------------------------------------------------------------------------------
#
#                                   KZFs Pipeline    
#
#   -----------------------------------------------------------------------------------------------------------

An pipeline designed to run RNA editing analysis on KZFs.

the pipeline input is sra accessions list or folder with fastq files

every .nf file have compatible .config file in the Configs folder

### run the nf file with:

nextflow -c file.config run file.nf 

set your deserved genome and dependencies by :
-profile
for example: -profile hg38

if you run parts that including STAR subpipeline you also should run 
----genome_length ${length}

we supporting length of 50,75,100,125,150

## for running the whole pipeline run:



# installation:

    install Nextflow and all its dependency https://nextflow.io/ 

    run init.sh file - it will download all the heavy resources from the web and will create salmon and star indexes
    !!!!  run init.sh wil take many hours !!!!


##    If you not run on levanon lab's servers:
###     you should have python 3.9 being able to run with the command "python3" and the following python modulus:
    pandas
    shutil
    subprocess
    gzip
    argparse
    re
    csv
    inspect
    json
    logging
    subprocess
    traceback
### you also need R being able to run with "Rscript" command, and the following libraries:
    argparse, 
    data.table, 
    dplyr, 
    dtplyr, 
    log4r, 
    purrr, 
