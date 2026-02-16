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

params.sample_table = false
params.targets_flanks = false
params.targets_regions = false
params.target_cpgs = false
params.input_type = false
params.genome = false
params.db = params.genome ? params.genomes[ params.genome ].db ?: false : false
params.multiqc_config = "${projectDir}/modules/multiqc/multiqc_config.yaml"
params.scripts = "${projectDir}/scripts/"

/*
================================================================================
Include modules to main pipeline
================================================================================
*/

/*
================================================================================
Include functions to main pipeline
================================================================================
*/

include { createInputChannel } from './functions/main.nf'

/*
================================================================================
Include subworkflows to main pipeline
================================================================================
*/

include { PREPROCESS } from './subworkflows/preprocess'
include { ALIGN } from './subworkflows/align'
include { READ_LEVEL_METHYLATION } from './subworkflows/read_level_methylation'
include { OLD_ASSAY_CALCULATION } from './subworkflows/old_assay_calculation'
include { RUN_MULTIQC } from './subworkflows/run_multiqc'

/*
================================================================================
Workflow declaration
================================================================================
*/

workflow {
    // Create an empty channel for multiqc input
    multiqc_ch = channel.empty()

    // Ingest sample table to create input channel
    input_ch = createInputChannel(params.sample_table, params.input_type)

    PREPROCESS(input_ch)

    ALIGN(
        PREPROCESS.out.fastp
    )

    READ_LEVEL_METHYLATION(
        ALIGN.out.clipped_reads
    )

    OLD_ASSAY_CALCULATION(
        ALIGN.out.cpg,
        ALIGN.out.noncpg,
        READ_LEVEL_METHYLATION.out.ref
    )

    multiqc_ch = multiqc_ch.mix(PREPROCESS.out.mqc.collect(), ALIGN.out.mqc.collect()) // add meta files to multiqc channel

    RUN_MULTIQC(
        multiqc_ch
    )
}