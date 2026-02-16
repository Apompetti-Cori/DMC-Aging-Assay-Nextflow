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
params.pubdir = "feature_counts"

/*
================================================================================
Module declaration
================================================================================
*/

process FEATURE_COUNTS {

    maxForks 4
    memory '8 GB'
    cpus 4
    conda "bioconda::subread"

    // Set batch name and sample id to tag
    tag { meta.batch == '' ? "${meta.id}" : "${meta.batch}_${meta.id}" }

    // Custom error handling logic
    errorStrategy {

    if (task.exitStatus == 255) {
            log.warn "WARNING: Skipping ${task.name} for '${meta.id}' no paired-end reads were detected."
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
    each path(target_regions)

    output:
    tuple val(meta), path("*.featureCounts.txt"), emit: counts
    tuple val(meta), path("*.featureCounts.txt.summary"), emit: meta_files

    script:
    """
    # Convert target regions from BED to SAF format (and convert to 1-based start)
    awk 'BEGIN{OFS="\\t"; print "GeneID\\tChr\\tStart\\tEnd\\tStrand"} {print \$4, \$1, \$2+1, \$3, "."}' ${target_regions} > target_regions.saf

    # Run featureCounts
    featureCounts --fracOverlapFeature 0.8 --countReadPairs -p -a target_regions.saf -F SAF -o ${meta.id}.featureCounts.txt ${reads[0]}
    """

}