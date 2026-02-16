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

include { PREPROCESS_READS } from '../../modules/preprocess_reads/main.nf'
include { FASTP } from '../../modules/fastp/main.nf'
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

workflow PREPROCESS {

    take:
        input_ch

    main:

        // Create an empty channel for multiqc input
        multiqc_ch = channel.empty()

        // Preprocess the sample table to change the files listed inside. Concatenates any multilane files.
        PREPROCESS_READS(input_ch)

        // Trim and filter reads 
        FASTP(PREPROCESS_READS.out)
        multiqc_ch = multiqc_ch.mix(FASTP.out.meta_files) // add fastp meta files to multiqc channel

    emit:
        fastp = FASTP.out.reads
        mqc = multiqc_ch.collect()
         
}