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


include { MULTIQC } from '../../modules/multiqc/main.nf'
include { RLM } from '../../modules/rlm/main.nf'
include { FILTER_RLM } from '../../modules/filter_rlm/main.nf'
include { GENERATE_CPG_DIST } from '../../modules/generate_cpg_dist/main.nf'
include { CREATE_REFERENCE } from '../../modules/create_reference/main.nf'
include { CALCULATE_DISTANCE } from '../../modules/calculate_distance/main.nf'

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

workflow READ_LEVEL_METHYLATION {

    take:
        reads_ch

    main:

    // Create an empty channel for multiqc input
    multiqc_ch = channel.empty()

    RLM(
        reads_ch,
        channel.fromPath( "${params.db}" )
    )

    FILTER_RLM(
        RLM.out.rlm,
        channel.fromPath( "${params.target_regions}" ),
        channel.fromPath( "${params.scripts}" )
    )

    GENERATE_CPG_DIST(
        FILTER_RLM.out.rlm,
        channel.fromPath( "${params.scripts}" )
    )

    ref_ch = GENERATE_CPG_DIST.out.dist
        .collect(flat: false)
        .flatMap()
        .filter{meta, _dist ->
            meta.condition == "reference"
        }
        .map{_meta, dist -> dist }
        .collect()

    CREATE_REFERENCE(
        ref_ch,
        channel.fromPath( "${params.scripts}" )
    )

    CALCULATE_DISTANCE(
        GENERATE_CPG_DIST.out.dist,
        CREATE_REFERENCE.out.refdist,
        channel.fromPath( "${params.scripts}" )
    )

    emit:
        ref = CREATE_REFERENCE.out.refdist
        mqc = multiqc_ch.collect()

}