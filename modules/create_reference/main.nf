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
params.pubdir = "ref_dist"

/*
================================================================================
Module declaration
================================================================================
*/

process CREATE_REFERENCE {

    maxForks 4
    memory '8 GB'
    cpus 4

    publishDir "${params.outdir}", mode: 'link', saveAs: { 
      filename ->
        "${params.pubdir}/${filename}"
    }

    input:
    path(dist), stageAs: 'dist_files/?/*'
    each path(scripts)

    output:
    path("reference_cpg_dist.bed"), emit: refdist

    script:

    """
    Rscript ${scripts}/create_reference/script.R \
      --input "${dist.join(',')}" \
      --output reference_cpg_dist.bed
    """
}