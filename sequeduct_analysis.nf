#!/usr/bin/env nextflow

params.sample_sheet = ''  // CSV of the sample ~ barcode relations. Columns: Sample,Barcode_dir (may have other)
params.references = ''   // dir of reference sequence Genbank files. Filenames must match 'Sample' column entries

params.fastq_dir = 'fastq'  // The directory that contains the barcode directories of FASTQ files
params.barcode_prefix = 'barcode'  // The prefix for the individual barcode directories as output by the sequencer

quality_cutoff = 9
max_len_fraction = 1.5
min_length = 500  // nucleotides
max_length = 12000  // temporary value


Channel
    .fromPath(params.sample_sheet)
    .splitCsv(header: true)
    .unique { row -> row['Barcode_dir'] }  // a barcode may be present multiple times due to plasmid pooling
    .map { row -> row.Barcode_dir }
    .set { barcodes_ch  }


process runNanoFilt {
    publishDir 'results/dir2_analysis/n1_fastq_filtered', mode: 'copy'

    input:
        val barcode from barcodes_ch

    output:
        // path fastq_file into fastq_filtered_ch
        tuple val(barcode), path(fastq_file) into fastq_filtered_ch
        // stdout result

    script:
        barcode_path = params.fastq_dir + '/' + barcode
        fastq_file = barcode + '.fastq'  // need for output
        fastqDir = file(barcode_path)

        fastqFiles = fastqDir.listFiles()  // multiple FASTQ in each barcode
        fastqFilePaths = []
        fastqFileString = fastqFiles.join(' ')  // need as one string for cat
        """
        cat $fastqFileString | NanoFilt -l $min_length --maxlength $max_length -q $quality_cutoff > $fastq_file
        """
}

process runNanoPlot {
    publishDir 'results/dir2_analysis/n2_nanoplots', mode: 'copy'

    input:
        tuple val(barcode), file(fastq_file) from fastq_filtered_ch
        // val fastq from fastq_filtered_ch
        // may need to convert fastq to string ?

    output:
        path barcode into nanoplots
        stdout result

    script:
        // barcode_path = params.fastq_dir + '/' + barcode
        // fastqDir = file(barcode_path)

        // fastqFiles = fastqDir.listFiles()  // multiple FASTQ in each barcode
        // fastqFilePaths = []
        // fastqFileString = fastqFiles.join(' ')  // need as one string for NanoPlot

        """
        NanoPlot --raw --fastq $fastq_file -o $barcode
        """
}


result.view()
