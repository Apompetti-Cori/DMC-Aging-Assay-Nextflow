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
params.pubdir = "bs_efficiency"

/*
================================================================================
Module declaration
================================================================================
*/

process BS_EFFICIENCY {

    maxForks 4
    memory '8 GB'
    cpus 4

    // Set batch name and sample id to tag
    tag { meta.batch == '' ? "${meta.id}" : "${meta.batch}_${meta.id}" }

    // Custom error handling logic
    errorStrategy {
        if (task.exitStatus == 2) {
            log.warn "WARNING: Skipping ${task.name} for '${meta.id}' because bismark cov file not found."
            return 'ignore'
        } else {
            return 'terminate'
        }
    }

    // Check batch and save output accordingly
    publishDir "${params.outdir}", mode: 'link', saveAs: { 
      filename ->
        return meta.batch == '' ? "${meta.id}/${params.pubdir}/${filename}" : "${meta.batch}/${meta.id}/${params.pubdir}/${filename}"
    }

    input:
    tuple val(meta), path(chg), path(chh)
    each path(scripts)

    output:
    tuple val(meta.id), val(meta.batch), path("*bs_efficiency.tsv"), emit: bseff, optional: true

    script:

    """
    Rscript ${scripts}/bs_efficiency/script.R ${chg} ${chh} ${meta.id}_bs_efficiency.tsv
    """
}