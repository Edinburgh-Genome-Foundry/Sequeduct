#!/usr/bin/env nextflow

params.sample_sheet = ''  // CSV of the sample ~ barcode relations. Columns: Sample,Barcode_dir (may have other)

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
        path barcode_fastq into fastq_filtered
        stdout result

    script:
        barcode_path = params.fastq_dir + '/' + barcode
        barcode_fastq = barcode + '.fastq'  // need for output
        fastqDir = file(barcode_path)

        fastqFiles = fastqDir.listFiles()  // multiple FASTQ in each barcode
        fastqFilePaths = []
        fastqFileString = fastqFiles.join(' ')  // need as one string for cat
        """
        cat $fastqFileString | NanoFilt -l $min_length --maxlength $max_length -q $quality_cutoff > $barcode_fastq
        """
}

result.view()
