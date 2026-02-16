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
params.pubdir = "bismark_align"

/*
================================================================================
Module declaration
================================================================================
*/

process BISMARK_ALIGN {

    maxForks 4
    memory '8 GB'
    cpus 4
    conda "bioconda::bismark bioconda::samtools"

    // Set batch name and sample id to tag
    tag { meta.batch == '' ? "${meta.id}" : "${meta.batch}_${meta.id}" }

    // Check batch and save output accordingly (Don't save trimmed fastq files)
    publishDir "${params.outdir}", mode: 'link', saveAs: { 
      filename ->
        return meta.batch == '' ? "${meta.id}/${params.pubdir}/${filename}" : "${meta.batch}/${meta.id}/${params.pubdir}/${filename}"
    }

    input:
    tuple val(meta), path(reads)
    each path(db)

    output:
    tuple val(meta), path("*.sorted.bam*"), emit: reads
    tuple val(meta), path("${meta.id}_{pe,se}.bam"), emit: unsorted_reads
    tuple val(meta), path("*{.idxstat,.flagstat,[SP]E_report}.txt"), emit: meta_files

    script:

    """
    bismark \
      --bowtie2 \
      -p ${task.cpus} \
      --un \
      --genome ${db} \
      --basename ${meta.id} \
      -1 ${reads[0]} \
      -2 ${reads[1]}
      samtools sort ${meta.id}_pe.bam -O bam -o ${meta.id}.sorted.bam
      samtools index ${meta.id}.sorted.bam

      samtools flagstat ${meta.id}.sorted.bam > ${meta.id}.flagstat.txt
      samtools idxstats ${meta.id}.sorted.bam > ${meta.id}.idxstat.txt
    """

}