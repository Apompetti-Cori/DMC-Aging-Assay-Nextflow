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

params.outdir = "./data/nfoutput"
params.pubdir = "calc_summary"

/*
================================================================================
Module declaration
================================================================================
*/

process CALC_SUMMARY {

    maxForks 4
    memory '8 GB'
    cpus 4

    // Set batch name and sample id to tag
    tag { batch == '' ? "${id}" : "${batch}_${id}" }

    // Custom error handling logic
    errorStrategy {
        if (task.exitStatus == 2) {
            log.warn "WARNING: Skipping ${task.name} for '${id}' because output file not found."
            return 'ignore'
        } else {
            return 'terminate'
        }
    }

    storeDir "${launchDir}/.nextflow/store/${batch}/${id}/${params.pubdir}/"

    // Check batch and save output accordingly
    publishDir "${params.outdir}", mode: 'link', saveAs: { 
      filename ->
        return batch == '' ? "${id}/${params.pubdir}/${filename}" : "${batch}/${id}/${params.pubdir}/${filename}"
    }

    input:
    tuple val(id), val(batch), path(bseff), path(afreq)
    each path(ref)
    each path(scripts)

    output:
    tuple val(id), val(batch), path("*_meth_jsd.tsv"), emit: meth_jsd, optional: true
    tuple val(id), val(batch), path("*.rds"), emit: rds, optional: true

    script:

    """
    Rscript ${scripts}/calc_summary/script.R ${ref} ${afreq} ${bseff} ${id}_meth_jsd.tsv
    """
}