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

workflow RUN_MULTIQC {

    take:
        multiqc_ch

    main:
        // Prepare multiqc input channel
        multiqc_ch = multiqc_ch
            .collect(flat: true)
            .flatMap { items ->
                // 'items' is now a flat list [meta1, file1, meta2, file2, ...]
                // .collate(2) groups it into a nested list [[meta1, file1], [meta2, file2], ...]
                items.collate(2).collect { meta, files ->
                    def batch = meta.batch

                    return tuple(batch, files)
                }
            }
            .groupTuple(by: 0) // Group by batch
            .map { batch, file_lists ->
                tuple(batch, file_lists.flatten()) // flatten file lists
            }

        //multiqc_ch.view()

        // Run multiqc
        MULTIQC(
            multiqc_ch,
            channel.fromPath( params.multiqc_config )
        )
}