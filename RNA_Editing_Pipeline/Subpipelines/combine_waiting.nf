nextflow.enable.dsl=2

/*
***********************************************************
this is just stupid workflow for alowing synchronize between processes
*********************************************************
*/


process do_nothing {
    input:
    val first
    val second
    output:
    val res
    script:
    res = true
    """
    """
}

workflow COMBINE_WAITING {
    take:
    waiting_one
    waiting_two
    main:
    both_finished=do_nothing(waiting_one,waiting_two)
    emit:
    both_finished
}