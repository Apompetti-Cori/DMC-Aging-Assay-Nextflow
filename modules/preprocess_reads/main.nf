#!/usr/bin/env nextflow

/*
================================================================================
Coriell Institute for Medical Research

Contributors:
Anthony Pompetti <apompetti@coriell.org>
================================================================================
*/

/*
================================================================================
Enable Nextflow DSL2
================================================================================
*/
nextflow.enable.dsl=2

/*
================================================================================
Module declaration
================================================================================
*/
process PREPROCESS_READS {
    
    maxForks 4

    // Set batch name and sample id to tag
    tag { meta.batch == '' ? "${meta.id}" : "${meta.batch}_${meta.id}" }

    input:
    tuple val(meta), path(reads_1), path(reads_2)

    output:
    tuple val(meta), path("{*_mL_SE,*_sL_SE,*_mL_PE_*,*_sL_PE_*}.fq.gz")

    script:
    if ( meta.multi_lane ){
        if (meta.single_end){ // single end reads 
            """
            cat ${reads_1} > ${meta.id}_mL_SE.fq.gz
            """
        }
        else {
            """
            cat ${reads_1} > ${meta.id}_mL_PE_1.fq.gz
            cat ${reads_2} > ${meta.id}_mL_PE_2.fq.gz
            """
        }
    }
    else {
        if (meta.single_end){
            """
            mv ${reads_1} ${meta.id}_sL_SE.fq.gz
            """
        }
        else {
            """
            mv ${reads_1} ${meta.id}_sL_PE_1.fq.gz
            mv ${reads_2} ${meta.id}_sL_PE_2.fq.gz
            """
        }
    }
}