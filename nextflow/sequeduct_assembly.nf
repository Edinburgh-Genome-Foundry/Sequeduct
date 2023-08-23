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

///////////////////////////////////////////////////////////////////////////////
// De novo assembly from reads

process runNanoFilt {
    publishDir 'results/dir4_assembly/n1_fastq_filtered', mode: 'copy', pattern: '*.fastq'  // need only the fastq

    input:
        tuple val(barcode), path(barcode_path), val(fastq_files), val(length)

    output:
        tuple val(barcode), path(barcode_path), path(fastq_file), val(length)

    script:
        fastq_file = barcode + '.fastq'  // need for output

        fastqFileString = fastq_files.join(' ')  // need as one string for cat

        """
        cat $fastqFileString | NanoFilt -l $params.min_length -q $params.quality_cutoff > $fastq_file
        """
}


process assembleOnly {
    publishDir 'results/dir4_assembly/n2_de_novo_assembly', mode: 'copy'
    
    input:
        tuple val(barcode), path(barcode_path), path(fastq_file), val(length)
    output:
        tuple val(barcode), path(assembly_dir)
    script:
        assembly_dir = barcode + '_assembly'
        genomsize_param = 'genomeSize=' + length + 'k'
        """
        canu -p $params.assembly_prefix -d $assembly_dir $genomsize_param -nanopore $fastq_file
        """
}


process trimAssemblyOnly {
    publishDir 'results/dir4_assembly/n3_assembly_trimmed', mode: 'copy', pattern: '*_denovo.fasta'

    input:
        tuple val(barcode), path(assembly_dir)
    output:
        tuple val(barcode), path(assembly_dir), path(trimmed_denovo)
    
    script:
        trimmed_denovo = barcode + '_denovo.fasta'
        """
        trim_assembly.py "$assembly_dir" "$params.assembly_prefix" "$params.canu_postfix" "$trimmed_denovo" "$barcode"
        """
}


workflow assemble_denovo {
    take: entries_assembly_ch
    main:
        runNanoFilt(entries_assembly_ch)
        assembleOnly(runNanoFilt.out)
        trimAssemblyOnly(assembleOnly.out)
}
