# DMC-Aging-Assay-Nextflow
Code availability for DMC Aging Assay

# How to run
- In order to run this pipeline, populate the ![sample table](sample_table_template.csv) with the fastq files you will be using.
    - First, fill in the sample_id (unique to a run), sample (unique across runs), and batch (relative to the sample_id).
    - Second, fill the r[0-9]_L[0-9] columns with the fastqs assigned to the sample_id
        - If single end simply fill in r1_L1
        - If paired end fill in r1_L1 and r2_L1
        - If paired end multi lane fill in r1_L[0-9] and r2_L[0-9]
    - Third, fill in condition of the sample
        - If marked as reference, the sample will be included in the reference distribution for calculating JSD
        - If marked as sample, the sample will be excluded from the reference
        - All samples regardless of condition get compared to the reference distribution

- Configure the pipeline to include the genome you will be aligning the fastq files to
    - Refer to ![nextflow.config](nextflow.config)
        - Fill in bismark index formatted how 'hg19' is
        - Label your genome and reference it in the ![nextflow shell script](nextflow_template.sh).

- Supply regions where you are targeting CpGs for JSD calculation in bed format.
    - also supply flanks of these regions in bed format if you wish to clip the ends of the reads (optional)

- Supply input file classification (pipeline has only been configured to run with fastq files)
    - Indicate "fastq" or "bam" ("fastq" is default)

- Supply the cpgs that will be used to calculate JSD in bed format
    - This will be used to ensure that all CpGs are covered by a read

- Run nextflow shell script to execute Nextflow Pipeline