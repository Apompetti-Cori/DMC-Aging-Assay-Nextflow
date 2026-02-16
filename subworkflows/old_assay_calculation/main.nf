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
include { BS_EFFICIENCY } from '../../modules/bs_efficiency/main.nf'
include { ALLELE_FREQUENCY } from '../../modules/allele_frequency/main.nf'
include { CALC_SUMMARY } from '../../modules/calc_summary/main.nf'

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

workflow OLD_ASSAY_CALCULATION {

    take:
        cpg
        noncpg
        ref

    main:

    // Calculate bisulfite efficiency
    def bseff_ch = BS_EFFICIENCY(
        noncpg,
        channel.fromPath( "${params.scripts}" )
    )

    // Calculate allele frequency
    def afreq_ch = ALLELE_FREQUENCY(
        cpg,
        channel.fromPath( "${params.target_cpgs}" ),
        channel.fromPath( "${params.scripts}" )
    )

    // Join afreq and bseff channels
    def ch_joined = bseff_ch
        .map { sample, batch, efficiency_file -> 
            tuple([sample, batch], efficiency_file)
        }
        .join(
            afreq_ch.map { sample, batch, allele_file ->
                tuple([sample, batch], allele_file)
            }
        )
        .map { key, efficiency_file, allele_file ->
            tuple(key[0], key[1], efficiency_file, allele_file)
        }

    //ch_joined.view()
    // Calculate summary statistics
    CALC_SUMMARY(
        ch_joined,
        ref,
        channel.fromPath( "${params.scripts}" )
    )

    // Create an empty channel for multiqc input
    multiqc_ch = channel.empty()

    emit:
        mqc = multiqc_ch.collect()

}