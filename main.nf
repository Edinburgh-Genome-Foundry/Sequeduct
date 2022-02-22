#!/usr/bin/env nextflow
nextflow.enable.dsl=2

params.sample_sheet = ''  // CSV of the sample ~ barcode relations. Columns: Sample,Barcode_dir (may have other)
params.fastq_dir = 'fastq'  // The directory that contains the barcode directories of FASTQ files
params.barcode_prefix = 'barcode'  // The prefix for the individual barcode directories as output by the sequencer

//
params.reference_dir = ''  // dir of reference sequence Genbank files. Filenames (without extension) must match 'Sample' column entries

params.max_len_fraction = 1.5  // For calculating max length of filtered FASTQ reads, to avoid plasmid dimers
params.projectname = 'Noname'  // display in PDF

params.quality_cutoff = 9  // For NanoFilt
params.min_length = 500  // nucleotides

//

include { preview_workflow } from "$projectDir/nextflow/sequeduct_preview.nf"
include { analysis_workflow } from "$projectDir/nextflow/sequeduct_analysis.nf"


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

workflow analysis {

    Channel
    .fromPath(params.sample_sheet)
    .splitCsv(header: true)
    .map { row -> 
        def barcode_dir = row['Barcode_dir']  // a barcode is present only once -- no pooling
        def sample = row['Sample']  // a sample may be present multiple times, in different barcodes
        def entry = "${barcode_dir}_${sample}"  // unique key for each sample sheet entry
        def genbank_path = file("${params.reference_dir}/${sample}.gb")
        def barcode_path = file("${params.fastq_dir}/${barcode_dir}")
        def fastq_files = barcode_path.listFiles()  // multiple FASTQ in each barcode

        return [entry, barcode_dir, barcode_path, fastq_files, sample, genbank_path]
        }
    .set { entries_ch }

    Channel
        .fromPath(params.sample_sheet)
        .splitCsv(header: true)
        .unique { row -> row['Sample'] }
        .map { row -> file(params.reference_dir + '/' + row['Sample'] + '.gb') }
        .set { genbank_ch }

    analysis_workflow(entries_ch, genbank_ch)
}