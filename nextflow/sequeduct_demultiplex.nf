#!/usr/bin/env nextflow
// Copyright 2025 Edinburgh Genome Foundry, University of Edinburgh
//
// This file is part of Sequeduct.
//
// Sequeduct is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
//
// Sequeduct is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License along with Sequeduct. If not, see <https:www.gnu.org/licenses/>.

nextflow.enable.dsl=2


process createFASTA {
    input:
        tuple val(barcode), file(barcode_path), val(fastq_files), val(sample), file(genbank_paths)
    output:
        tuple val(barcode), file(barcode_path), val(fastq_files), val(sample), path(barcode_fasta)
    script:
        barcode_fasta = barcode + '.fa'
        gbFileString = genbank_paths.join(' ')
        """
        create_2x_fasta.py "$barcode_fasta" $gbFileString
        """
}

process alignMultiplexReads {
    input:
        tuple val(barcode), file(barcode_path), val(fastq_files), val(sample), path(barcode_fasta)
    output:
        tuple val(barcode), file(barcode_path), val(fastq_files), val(sample), path(paf_file)
    script:
        fastqFileString = fastq_files.join(' ')
        sam_file = barcode + '.sam'
        paf_file = barcode + '.paf'
        """
        cat $fastqFileString | \
        minimap2 --secondary=no -ax map-ont $barcode_fasta - > $sam_file
        paftools.js sam2paf $sam_file > $paf_file
        """
}

process createSubDirs {
   publishDir 'results/dir0_demultiplex/', mode: 'copy'
   input:
        tuple val(barcode), file(barcode_path), val(fastq_files), val(sample), path(paf_file)
    output:
        path "split_fastq/*"
    script:
        fastqFileString = fastq_files.join(' ')
        split_fastq_dir = "split_fastq"
        """
        demultiplex_data.py $split_fastq_dir $fastqFileString
        """
}

// process createSampleSheet {

// }

workflow demultiplex_workflow {
    take:
        multiplex_ch

    main:
        createFASTA(multiplex_ch)
        alignMultiplexReads(createFASTA.out)
        createSubDirs(alignMultiplexReads.out)
}
