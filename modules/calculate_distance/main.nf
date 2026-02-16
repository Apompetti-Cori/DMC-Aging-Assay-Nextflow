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
params.pubdir = "calc_dist"

/*
================================================================================
Module declaration
================================================================================
*/

process CALCULATE_DISTANCE {

    maxForks 4
    memory '8 GB'
    cpus 4

    // Check batch and save output accordingly
    publishDir "${params.outdir}", mode: 'link', saveAs: { 
      filename ->
        return meta.batch == '' ? "${meta.id}/${params.pubdir}/${filename}" : "${meta.batch}/${meta.id}/${params.pubdir}/${filename}"
    }

    input:
    tuple val(meta), path(dist)
    each path(refdist)
    each path(scripts)

    output:
    path("*_distance.bed"), emit: refdist

    script:

    """
    Rscript ${scripts}/calculate_distance/script.R \
      --input "${dist}" \
      --reference "${refdist}" \
      --output "${meta.sample}_distance.bed"
    """
}