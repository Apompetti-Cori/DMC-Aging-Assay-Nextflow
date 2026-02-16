#!/bin/bash

nextflow \
-log "./data/.pipeline_info/nextflow.log" \
run "path to main.nf file describing the pipeline" \
-with-dag "./data/.pipeline_info/pipeline_dag.svg" \
-with-report "./data/.pipeline_info/execution_report.html" \
-with-trace "./data/.pipeline_info/execution_trace.txt" \
-resume \
--sample_table "path to sample table that will be used for input fastq files" \
--input_type "fastq" \
--target_flanks "flanks of the target regions you want to use to trim reads" \
--target_regions "regions of the targets in amplicon assay" \
--target_cpgs "coordinates of the cpgs used in amplicon assay" \
--genome "genome you wish to align to"