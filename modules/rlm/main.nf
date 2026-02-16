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
Configurable variables for module
================================================================================
*/

/*
================================================================================
Module declaration
================================================================================
*/

process RLM {

    maxForks 4
    memory '8 GB'
    cpus 4

    // Set batch name and sample id to tag
    tag { meta.batch == '' ? "${meta.id}" : "${meta.batch}_${meta.id}" }

    input:
    tuple val(meta), path(reads)
    each path(db)

    output:
    tuple val(meta), path("*_single_read.bed"), emit: rlm

    script:

    """
    RLM -b ${reads[0]} \
      -r ${db}/hg19lambda.fa \
      --mode PE -a bismark \
      -c 1 \
      -q 10 \
      -s single_read \
      -o ${meta.id}_single_read.bed
    """

}