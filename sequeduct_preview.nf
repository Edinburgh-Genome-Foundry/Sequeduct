#!/usr/bin/env nextflow

params.sample_sheet = ''  // CSV of the sample ~ barcode relations. Columns: Sample,Barcode_dir (may have other)
params.fastq_dir = 'fastq'  // The directory that contains the barcode directories of FASTQ files
params.barcode_prefix = 'barcode'  // The prefix for the individual barcode directories as output by the sequencer

Channel
    .fromPath(params.sample_sheet)
    .splitCsv(header: true)
    .unique { row -> row['Barcode_dir'] }  // a barcode may be present multiple times due to plasmid pooling
    .map { row -> row.Barcode_dir }
    .set { barcodes_ch }

process runNanoPlot {
    publishDir 'results/dir1_preview', mode: 'copy'

    input:
        val barcode from barcodes_ch

    output:
        path barcode into nanoplots
        stdout result

    script:
        barcode_path = params.fastq_dir + '/' + barcode
        fastqDir = file(barcode_path)

        fastqFiles = fastqDir.listFiles()  // multiple FASTQ in each barcode
        fastqFilePaths = []
        fastqFileString = fastqFiles.join(' ')  // need as one string for NanoPlot

        """
        NanoPlot --raw --fastq $fastqFileString -o $barcode
        """
}

result.view()
