# DMC-Aging-Assay-Nextflow
Code availability for DMC Aging Assay

# How to run
- In order to run this pipeline, populate the sample table with the fastq files you will be using.
    - First, fill in the sample_id (unique to a run), sample (unique across runs), and batch (relative to the sample_id).
    - Second, fill the r[0-9]_L[0-9] columns with the fastqs assigned to the sample_id
        - If single end simply fill in r1_L1
        - If paired end fill in r1_L1 and r2_L1
        - If paired end multi lane fill in r1_L[0-9] and r2_L[0-9]
    - Third, fill in condition of the sample
        - If marked as reference, the sample will be included in the reference distribution for calculating JSD
        - If marked as sample, the sample will be excluded from the reference
        - All samples regardless of condition get compared to the reference distribution