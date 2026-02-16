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
Function declaration
================================================================================
*/
def createInputChannel(String sample_table, String input_type) {
    // Channel for the samplesheet
    ch_samplesheet = Channel.fromPath(sample_table)

    if(input_type == "fastq"){
        // Parse it line by line
        input_ch = ch_samplesheet.splitCsv(header:true).filter{
            // Filter out observations without files
            it['r1_L1'] != ''
        }.map{
            // This is the read1 and read2 entry
            r1_L1 = it['r1_L1']
            r1_L2 = it['r1_L2']
            r1_L3 = it['r1_L3']
            r1_L4 = it['r1_L4']
            r2_L1 = it['r2_L1']
            r2_L2 = it['r2_L2']
            r2_L3 = it['r2_L3']
            r2_L4 = it['r2_L4']
            merge_id = it['merge_id']
            merge_batch = it['merge_batch']
            condition = it['condition']

            // Detect wiether single-end or paired-end
            is_singleEnd = r2_L1.toString()=='' ? true : false
            is_multiLane = r1_L2.toString()=='' ? false : true
            is_emptyObsv = r1_L1.toString()=='' ? true : false
            merge = merge_id.toString()=='' ? false : true
            
            // The "meta" map, which is a Nextflow/Groovy map with id (the sample name) and a single_end logical entry
            meta = [id: it['sample_id'], sample: it['sample'], batch: it['batch'], merge_id: merge_id, merge_batch: merge_batch, merge: merge, single_end: is_singleEnd, multi_lane: is_multiLane, empty_obsv: is_emptyObsv, condition: condition]
            
            // We return a nested map, the first entry is the meta map, the second one is the read(s)
            if ( is_singleEnd ){
                if ( r1_L4.toString()!='' ){
                    [meta, [r1_L1, r1_L2, r1_L3, r1_L4], []]
                }
                else if ( r1_L3.toString()!='' ){
                    [meta, [r1_L1, r1_L2, r1_L3], []]
                }
                else if ( r1_L2.toString()!='' ){
                    [meta, [r1_L1, r1_L2], []]
                }
                else {
                    [meta, [r1_L1], []]
                }
            }
            else{
                if (r2_L4.toString()!=''){
                    [meta, [r1_L1, r1_L2, r1_L3, r1_L4], [r2_L1, r2_L2, r2_L3, r2_L4]]
                }
                else if ( r2_L3.toString()!='' ){
                    [meta, [r1_L1, r1_L2, r1_L3], [r2_L1, r2_L2, r2_L3]]
                }
                else if ( r2_L2.toString()!='' ){
                    [meta, [r1_L1, r1_L2], [r2_L1, r2_L2]]
                }
                else {
                    [meta, [r1_L1], [r2_L1]]
                }
            }
        }
    }
    if(input_type == "bam"){
        // Parse it line by line
        input_ch = ch_samplesheet.splitCsv(header:true).filter{
            // Filter out observations without files
            it['r1_L1'] != ''
        }.map{
        
        // Look at the column r1_L1 and use these files as input
        input = it['r1_L1']
        
        // The "meta" map, which is a Nextflow/Groovy map with id (the sample name) and other metadata
        meta = [id: it['sample_name'], batch: it['batch']]
        
        [meta, input]
        }
    }

    return input_ch
}

def groupSamplesByMeta(input_ch){

    // Group samples by the "sample" entry in the meta map
    def grouped_ch = input_ch
        .flatMap()
        .map{ meta, reads -> return [meta.sample, meta.batch, ["id": meta.id, "condition": meta.condition, "reads": reads]] }
        .groupTuple(by: [0,1])
        .map{ id, batch, entries_list ->
            // Initialize empty lists for input and cut reads for the current sample
            def input_reads
            def cut_reads
            def input_id
            def cut_id

            // Iterate through the list of entries for this sample
            entries_list.each { entry ->
                // Check the condition and add to the appropriate list
                if (entry.condition == 'input') {
                    input_reads = entry.reads
                    input_id = entry.id
                }

                if (entry.condition == 'cut') {
                    cut_reads = entry.reads
                    cut_id = entry.id
                }
            }

        return [["id": id, "batch": batch], input_reads, cut_reads, input_id, cut_id]
        }

    return grouped_ch
}

def collectCoverage(input_ch){

    // Collect all BAM files from the input channel
    def bw_files = input_ch
        .flatMap()
        .map{ meta, reads ->
        return [meta.batch, meta.condition, meta.sample, reads] }
        .groupTuple(by: 0)

    return bw_files
}