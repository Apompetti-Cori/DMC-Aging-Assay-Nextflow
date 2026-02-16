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

params.outdir = "./nfoutput"
params.pubdir = "multiqc"

/*
================================================================================
Module declaration
================================================================================
*/


process MULTIQC {

    memory '8 GB'
    cpus 1

    // Set batch name and sample id to tag
    tag { batch == '' ?: "${batch}" }

    conda 'bioconda::multiqc=1.25.1'
    
    // Check batch and save output accordingly
    publishDir "${params.outdir}", mode: 'link', saveAs: {
      filename ->
        return batch == '' ? "${params.pubdir}/${filename}" : "${params.pubdir}/${batch}/${filename}"
    }

    input:
    tuple val(batch), path('multiqc_input/?/*')
    each path(multiqc_config)

    output:
    path("*.html"), emit: html
    path("*_data"), emit: data, type: 'dir'

    script:
    if( batch != '' ) {
        """
        multiqc multiqc_input/ \
            --fn_as_s_name \
            --config ${multiqc_config} \
            --title "MultiQC Report ${batch}" \
            --filename ${batch}.multiqc_report.html
        """
    } else {
        """
        multiqc multiqc_input/ \
            --fn_as_s_name \
            --config ${multiqc_config} \
            --title "MultiQC Report" \
            --filename multiqc_report.html
        """
    }

}