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


///////////////////////////////////////////////////////////////////////////////

// Multiplex workflow processes:
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
    publishDir 'results/dir0_demultiplex', mode: 'copy', pattern: "$split_fastq_dir/*"
    input:
        tuple val(barcode), file(barcode_path), val(fastq_files), val(sample), path(paf_file)
    output:
        path "${split_fastq_dir}/*"
        path samplesheet_entries_csv, emit: samplesheet_entries_csv_ch
    script:
        fastqFileString = fastq_files.join(' ')
        split_fastq_dir = "fastq_split"
        samplesheet_entries_csv = "samplesheet_entries.csv"  // for use in the analysis pipeline
        """
        demultiplex_data.py $split_fastq_dir $barcode $paf_file $fastqFileString
        """
}

process createMultiPlexSampleSheet {
    // publishDir 'results/dir0_demultiplex', mode: 'copy'
    input:
        path(samplesheet_entries_csv)
    output:
        path(samplesheet_split_csv)
    script:
        samplesheet_split_csv = "samplesheet_split.csv"
        """
        (echo "Sample,Barcode_dir"; cat "$samplesheet_entries_csv") > $samplesheet_split_csv
        """
}

///////////////////////////////////////////////////////////////////////////////

// Singleplex workflow processes:
process copyFastqDir {
    publishDir 'results/dir0_demultiplex/fastq_split', mode: 'copy', pattern: "${barcode}", enabled: params.singleplex_out
    input:
        tuple val(barcode), file(barcode_path), val(fastq_files), val(sample), file(genbank_paths)
    output:
        file(barcode_path)
        path(samplesheet_csv), emit: singleplex_samplesheet_ch
    script:
        samplesheet_csv = "singleplex_samplesheet.csv"
        """
        echo "$sample,$barcode" >> $samplesheet_csv
        """
}

process createSinglePlexSampleSheet {
    // publishDir 'results/dir0_demultiplex', mode: 'copy'
    input:
        path(singleplex_samplesheet_ch)
    output:
        path(singleplex_samplesheet_ch)
    script:
        """
        """
}


///////////////////////////////////////////////////////////////////////////////

// Samplesheet workflow

process combineSampleSheets {
    publishDir 'results/dir0_demultiplex', mode: 'copy'
    input:
        path(singleplex_samplesheet_ch)
        path(multiplex_samplesheet_ch)
    output:
        path(combined_samplesheet_ch)
    script:
        combined_samplesheet_ch = "samplesheet.csv"
        """
        cat $multiplex_samplesheet_ch $singleplex_samplesheet_ch > $combined_samplesheet_ch
        """
}


///////////////////////////////////////////////////////////////////////////////

// Workflows

workflow demultiplex_workflow {
    take:
        multiplex_ch
    main:
        createFASTA(multiplex_ch)
        alignMultiplexReads(createFASTA.out)
        createSubDirs(alignMultiplexReads.out)
        createMultiPlexSampleSheet(createSubDirs.out.samplesheet_entries_csv_ch.collectFile())
    emit:
        createMultiPlexSampleSheet.out
}

workflow singleplex_workflow {
    take:
        singleplex_ch
    main:
        copyFastqDir(singleplex_ch)
        createSinglePlexSampleSheet(copyFastqDir.out.singleplex_samplesheet_ch.collectFile())
    emit:
        createSinglePlexSampleSheet.out
}

workflow combine_samplesheets {
    take:
        singleplexsamplesheet_ch
        multiplexsamplesheet_ch
    main:
        combineSampleSheets(singleplexsamplesheet_ch, multiplexsamplesheet_ch)
}
