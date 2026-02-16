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

process CLIP_READS {

    maxForks 4
    memory '8 GB'
    cpus 4
    conda "bioconda::bismark bioconda::samtools bioconda::bedtools"

    // Set batch name and sample id to tag
    tag { meta.batch == '' ? "${meta.id}" : "${meta.batch}_${meta.id}" }

    input:
    tuple val(meta), path(reads)
    each path(target_regions)
    each path(target_flanks)
    each path(db)

    output:
    tuple val(meta), path("*.clipped.sorted.bam*"), emit: reads

    script:

    """
    samtools ampliconclip \
      -b ${target_flanks} \
      -o ${meta.id}.clipped.bam \
      --hard-clip \
      --both-ends \
      ${reads[0]}

    samtools sort -o ${meta.id}.clipped.sorted.bam ${meta.id}.clipped.bam
    samtools index ${meta.id}.clipped.sorted.bam
    """

}