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

params.outdir = "data/nfoutput"
params.pubdir = "fastp"

/*
================================================================================
Module declaration
================================================================================
*/

process FASTP {

    maxForks 4
    memory '8 GB'
    cpus 4

    // Set batch name and sample id to tag
    tag { meta.batch == '' ? "${meta.id}" : "${meta.batch}_${meta.id}" }

    // Check batch and save output accordingly (Don't save trimmed fastq files)
    publishDir "${params.outdir}", mode: 'link', saveAs: { 
      filename ->
        if (filename.endsWith(".gz")) {
            return null
        } else {
            return meta.batch == '' ? "${meta.id}/${params.pubdir}/${filename}" : "${meta.batch}/${meta.id}/${params.pubdir}/${filename}"
        }
    }

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("*.fq.gz"), emit: reads
    tuple val(meta), path("*.{html,json}"), emit: meta_files

    script:

    if ( meta.single_end ){

        """
        fastp \
        --in1 ${reads} \
        --out1 ${meta.id}.trimmed.fq.gz \
        -h ${meta.id}.fastp.html \
        -j ${meta.id}.fastp.json \
        -w ${task.cpus}
        """

    }
    else {

        """
        fastp \
        --in1 ${reads[0]} \
        --in2 ${reads[1]} \
        --out1 ${meta.id}.trimmed.1.fq.gz \
        --out2 ${meta.id}.trimmed.2.fq.gz \
        -h ${meta.id}.fastp.html \
        -j ${meta.id}.fastp.json \
        -w ${task.cpus}
        """

    }

}