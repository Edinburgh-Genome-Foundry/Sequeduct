#!/usr/bin/env nextflow
nextflow.enable.dsl=2

process runNanoPlot {
    publishDir 'results/dir1_preview', mode: 'copy'

    input:
        tuple val(barcode), file(barcode_path), val(fastq_files)
    output:
        path barcode

    script:
        fastq_file_string = fastq_files.join(' ')  // need as one string for cat
        """
        NanoPlot --raw --fastq $fastq_file_string -o $barcode
        """
}

workflow preview_workflow {
    take: input_ch
    main:
        runNanoPlot(input_ch)
}

