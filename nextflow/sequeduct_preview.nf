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

