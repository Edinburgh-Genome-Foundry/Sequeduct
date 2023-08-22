#!/usr/bin/env nextflow
// Copyright 2021 Edinburgh Genome Foundry, University of Edinburgh
//
// This file is part of Sequeduct.
//
// Sequeduct is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
//
// Sequeduct is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License along with Sequeduct. If not, see <https:www.gnu.org/licenses/>.

nextflow.enable.dsl=2


include { preview_workflow } from "$projectDir/nextflow/sequeduct_preview.nf"
include { analysis_workflow } from "$projectDir/nextflow/sequeduct_analysis.nf"
include { review_consensus } from "$projectDir/nextflow/sequeduct_review.nf"
include { review_denovo } from "$projectDir/nextflow/sequeduct_review.nf"
include { assemble_denovo } from "$projectDir/nextflow/sequeduct_review.nf"


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

workflow review {

    Channel
        .fromPath(params.results_csv)
        .splitCsv(header: true)
        .filter { item -> item[params.consensus_columname] == params.consensus_true}
        
        .map { row -> 
            def barcode = row['Barcode']  // a barcode is present only once -- no pooling
            def sample = row['Sample']  // a sample may be present multiple times, in different barcodes
            def entry = "${barcode}_${sample}"  // unique key for each sample sheet entry
            def result = row["Result"]
            def genbank_path = file("${params.reference_dir}/${sample}.gb")
            // Note: sample name matches filename with ".gb" extension
            def consensus_path = file("${params.consensus_dir}/${entry}_consensus.fa")
            return [entry, barcode, sample, result, genbank_path, consensus_path]
            }
        .set { entries_ch }

    Channel
        .fromPath(params.results_csv)
        .splitCsv(header: true)
        .filter { item -> item[params.denovo_columname] == params.denovo_true}
        
        .map { row -> 
            def barcode = row['Barcode']  // a barcode is present only once -- no pooling
            def sample = row['Sample']  // a sample may be present multiple times, in different barcodes
            def entry = "${barcode}_${sample}"  // unique key for each sample sheet entry
            def result = row["Result"]
            def genbank_path = file("${params.reference_dir}/${sample}.gb")
            def fastq_path = file("${params.fastq_filtered_dir}/${barcode}.fastq")
            return [entry, barcode, sample, result, genbank_path, fastq_path]
            }
        .set { entries_de_novo_ch }


    review_consensus(entries_ch)

    review_denovo(entries_de_novo_ch)
}


workflow assembly {
    Channel
        .fromPath(params.results_csv)
        .splitCsv(header: true)
        
        .map { row -> 
            def barcode = row['Barcode']  // a barcode is present only once -- no pooling
            def length = row['Length']  // estimated length in kbp for canu genomeSize parameter
            def fastq_path = file("${params.fastq_filtered_dir}/${barcode}.fastq")
            return [barcode, fastq_path, length]
            }
        .set { entries_assembly_ch }

   assemble_denovo(entries_assembly_ch)
}