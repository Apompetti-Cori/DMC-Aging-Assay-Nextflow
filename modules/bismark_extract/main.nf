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
params.pubdir = "bismark_extract"

/*
================================================================================
Module declaration
================================================================================
*/

process BISMARK_EXTRACT {

    maxForks 4
    memory '8 GB'
    cpus 4
    conda "bioconda::bismark bioconda::samtools"

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
    tuple val(meta), path(reads)
    each path(db)

    output:
    tuple val(meta), path("*.bismark.cov.gz"), emit: cov, optional: true
    tuple val(meta), path("*.bedGraph.gz"), emit: bedgraph, optional: true
    tuple val(meta), path("*{_splitting_report,M-bias}.txt"), emit: meta_files
    tuple val(meta), path("CpG_context*.txt.gz"), emit: cpg
    tuple val(meta), path("CHG_context*.txt.gz"), path("CHH_context*.txt.gz"), emit: noncpg

    script:

    """
    bismark_methylation_extractor \
      --multicore ${task.cpus} \
      --bedGraph \
      --counts \
      --gzip \
      --comprehensive \
      ${reads}

    files=(*.bismark.cov.gz)
    if [ -e "\${files[0]}" ]; then
        echo "Found .bismark.cov.gz files"
    else
        echo "No .bismark.cov.gz files found"
        exit 2
    fi
    """
}