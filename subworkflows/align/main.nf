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
Configurable variables for pipeline
================================================================================
*/
params.multiqc_config = "${projectDir}/modules/multiqc/multiqc_config.yaml"

/*
================================================================================
Include modules to main pipeline
================================================================================
*/

include { BISMARK_ALIGN } from '../../modules/bismark_align/main.nf'
include { FEATURE_COUNTS } from '../../modules/feature_counts/main.nf'
include { BISMARK_EXTRACT } from '../../modules/bismark_extract/main.nf'
include { CLIP_READS } from '../../modules/clip_reads/main.nf'
include { MULTIQC } from '../../modules/multiqc/main.nf'

/*
================================================================================
Include functions to main pipeline
================================================================================
*/

/*
================================================================================
Workflow declaration
================================================================================
*/

workflow ALIGN {

    take:
        reads_ch

    main:

        // Create an empty channel for multiqc input
        multiqc_ch = channel.empty()

        // Align trimmed reads to reference using bwa-mem2
        BISMARK_ALIGN(
            reads_ch,
            channel.fromPath( "${params.db}" )
        )
        multiqc_ch = multiqc_ch.mix(BISMARK_ALIGN.out.meta_files) // add samtools meta files to multiqc channel

        // Count number of reads that aligned to targeted regions
        FEATURE_COUNTS(
            BISMARK_ALIGN.out.reads,
            channel.fromPath( "${params.target_regions}" )
        )
        multiqc_ch = multiqc_ch.mix(FEATURE_COUNTS.out.meta_files) // add featureCounts meta files to multiqc channel

        // Trim reads to keep only sections desired
        CLIP_READS(
            BISMARK_ALIGN.out.reads,
            channel.fromPath( "${params.target_regions}" ),
            channel.fromPath( "${params.target_flanks}" ),
            channel.fromPath( "${params.db}" )
        )

        // Extract methylation calls from clipped reads
        BISMARK_EXTRACT(
            BISMARK_ALIGN.out.unsorted_reads,
            channel.fromPath( "${params.db}" )
        )
        multiqc_ch = multiqc_ch.mix(BISMARK_EXTRACT.out.meta_files) // add bismark meta files to multiqc channel

    emit:
        clipped_reads = CLIP_READS.out.reads
        cpg = BISMARK_EXTRACT.out.cpg
        noncpg = BISMARK_EXTRACT.out.noncpg
        mqc = multiqc_ch.collect()

}