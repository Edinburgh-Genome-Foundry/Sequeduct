#!/usr/bin/env nextflow

params.fastq_dir = 'fastq/barcode01'

fastqDir = file(params.fastq_dir)
fastqFiles = fastqDir.list()
fastqFilePaths = []
for( i in fastqFiles ) {
    fastqFilePaths += PWD + "/" + params.fastq_dir + "/" + i // join path
}
fastqFileString = fastqFilePaths.join(' ')

process runNanoPlot {

	publishDir 'results/', mode: 'copy'

    // input:

    output:
    path 'nanoplots' into nanoplots
    stdout result

    """
    NanoPlot --raw --fastq $fastqFileString -o nanoplots
    # echo "$fastqFileString"
    """
}

result.view()
