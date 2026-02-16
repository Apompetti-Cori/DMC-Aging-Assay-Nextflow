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

process FILTER_RLM {

    maxForks 4
    memory '8 GB'
    cpus 4

    // Set batch name and sample id to tag
    tag { meta.batch == '' ? "${meta.id}" : "${meta.batch}_${meta.id}" }

    // Custom error handling logic
    errorStrategy {

        if (task.exitStatus == 42) {
            log.warn "WARNING: Skipping ${task.name} for '${meta.id}' because the input RLM file was empty."
            return 'ignore'
        } else if (task.exitStatus == 43) {
            log.warn "WARNING: Skipping ${task.name} for '${meta.id}' because the output RLM file was empty."
            return 'ignore'
        } else {
            return 'terminate'
        }
    }

    input:
    tuple val(meta), path(rlm)
    each path(targets)
    each path(scripts)

    output:
    tuple val(meta), path("*_single_read.filtered.bed"), emit: rlm

    script:

    """
    Rscript ${scripts}/filter_rlm/script.R \
      --bed_file ${targets} \
      --input ${rlm} \
      --output ${meta.id}_single_read.filtered.bed \
      --unannotated_output ${meta.id}_single_read.unannotated.bed
    """

}