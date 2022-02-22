#!/usr/bin/env nextflow
nextflow.enable.dsl=2

params.sample_sheet = ''  // CSV of the sample ~ barcode relations. Columns: Sample,Barcode_dir (may have other)
params.fastq_dir = 'fastq'  // The directory that contains the barcode directories of FASTQ files
params.barcode_prefix = 'barcode'  // The prefix for the individual barcode directories as output by the sequencer

include { preview_workflow } from "$projectDir/nextflow/sequeduct_preview.nf"


workflow preview {
    Channel
        .fromPath(params.sample_sheet)
        .splitCsv(header: true)
        .unique { row -> row['Barcode_dir'] }  // a barcode may be present multiple times due to plasmid pooling
        .map { row ->
            def barcode_out = row['Barcode_dir'] + '_plots'
            def barcode_dir = row['Barcode_dir']  // a barcode is present only once -- no pooling
            def barcode_path = file("${params.fastq_dir}/${barcode_dir}")
            def fastq_files = barcode_path.listFiles()  // multiple FASTQ in each barcode
            return [barcode_out, barcode_path, fastq_files]
        }
        .set { input_ch }
    preview_workflow(input_ch)
}
