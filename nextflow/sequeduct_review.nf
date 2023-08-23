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
// Variant call consensus sequence review

process convertGenbank {

    input:
        tuple val(entry), val(barcode), val(sample), val(result), file(genbank_path), file(consensus_path)

    output:
        tuple val(entry), val(barcode), val(sample), val(result), file(genbank_path), path(sample_fasta), file(consensus_path)

    script:
        sample_fasta = sample + '.fa'

        """
        convert_genbank.py "$genbank_path" "$sample" "$sample_fasta" "$params.max_len_fraction" "none"
        """
}

process alignParts {
    publishDir 'results/dir3_review/n1_consensus_alignment', mode: 'copy', pattern: '*.paf'

    input:
        tuple val(entry), val(barcode), val(sample), val(result), val(genbank_path), path(sample_fasta), file(consensus_path)
        file parts_path

    output:
        tuple val(entry), val(barcode), val(sample), val(result), val(genbank_path), path(sample_fasta), file(consensus_path), path(paf), emit: alignment_ch

    script:
        paf = entry + '.paf'
        """
        cat $sample_fasta $parts_path | \
        minimap2 -cx asm5 $consensus_path - > $paf
        """
}

process writeCSV {
    input:
        tuple val(entry), val(barcode), val(sample), val(result), val(genbank_path), path(sample_fasta), path(consensus_path), path(paf)
    output:
        path samplesheet_csv, emit: samplesheet_csv_ch
        path paf, emit: paf_file_ch
        path consensus_path, emit: consensus_path_ch
    script:
        samplesheet_csv = "entries.csv"
        // order is important, see Python script:
        """
        echo "$params.projectname,$entry,$barcode,$sample,$result,$genbank_path,$sample_fasta,$consensus_path,$paf" >> $samplesheet_csv
        """    
}

process runReview {
    publishDir 'results/dir3_review/n2_consensus_results', mode: 'copy'

    input:
        file paf
        file consensus
        path samplesheet_csv
    output:
        tuple path(pdf_file), path(samplesheet_csv)
    script:
        pdf_file = "consensus_review.pdf"
        """
        review.py "$samplesheet_csv" "$params.plan_path" "$params.projectname" "$pdf_file"
        """
}


workflow review_consensus {
    take: entries_ch
    main:
        params.parts_path = file(params.all_parts)
        if( params.assembly_plan == "noplan" ) {
            params.plan_path = "noplan"
        }
        else {
            params.plan_path = file(params.assembly_plan)
        }
        convertGenbank(entries_ch)
        alignParts(convertGenbank.out, params.parts_path)
        writeCSV(alignParts.out.alignment_ch)
        runReview(writeCSV.out.paf_file_ch.collect(), writeCSV.out.consensus_path_ch.collect(), writeCSV.out.samplesheet_csv_ch.collectFile())
}

///////////////////////////////////////////////////////////////////////////////
// De novo assembly sequence review

process convertGenbank_de_novo {
    input:
        tuple val(entry), val(barcode), val(sample), val(result), file(genbank_path), file(fastq_path)

    output:
        tuple val(entry), val(barcode), val(sample), val(result), file(genbank_path), path(sample_fasta), stdout, path(fastq_path)  // stdout for seq length

    script:
        sample_fasta = sample + '.fa'

        """
        convert_genbank.py "$genbank_path" "$sample" "$sample_fasta" "$params.max_len_fraction" "canu"
        """
}

process assembleDeNovo {
    publishDir 'results/dir3_review/n3_de_novo_assembly', mode: 'copy'
    
    input:
        tuple val(entry), val(barcode), val(sample), val(result), file(genbank_path), path(sample_fasta), val(seq_length), path(fastq_path)
    output:
        tuple val(entry), val(barcode), val(sample), val(result), file(genbank_path), path(sample_fasta), path(assembly_dir)
    script:
        assembly_dir = barcode + '_assembly'
        genomsize_param = 'genomeSize=' + seq_length + 'k'
        """
        canu -p $params.assembly_prefix -d $assembly_dir $genomsize_param -nanopore $fastq_path
        """
}

process trimAssembly {
    publishDir 'results/dir3_review/n3_de_novo_assembly/trimmed', mode: 'copy', pattern: '*_denovo.fasta'

    input:
        tuple val(entry), val(barcode), val(sample), val(result), val(genbank_path), path(sample_fasta), val(assembly_dir)
    output:
        tuple val(entry), val(barcode), val(sample), val(result), val(genbank_path), path(sample_fasta), val(assembly_dir), path(trimmed_denovo)
    
    script:
        trimmed_denovo = barcode + '_denovo.fasta'
        """
        trim_assembly.py "$assembly_dir" "$params.assembly_prefix" "$params.canu_postfix" "$trimmed_denovo" "$barcode"
        """
}

process alignParts_de_novo {
    publishDir 'results/dir3_review/n4_de_novo_alignment', mode: 'copy', pattern: '*.paf'

    input:
        tuple val(entry), val(barcode), val(sample), val(result), val(genbank_path), path(sample_fasta), val(assembly_dir), path(trimmed_denovo)
        file parts_path

    output:
        tuple val(entry), val(barcode), val(sample), val(result), val(genbank_path), path(sample_fasta), val(assembly_dir), path(trimmed_denovo), path(paf)

    script:
        paf = entry + '.paf'
        """
        cat $sample_fasta $parts_path | \
        minimap2 -cx asm5 $trimmed_denovo - > $paf
        """
}

process writeCSV_de_novo {
    input:
        tuple val(entry), val(barcode), val(sample), val(result), file(genbank_path), path(sample_fasta), val(assembly_dir), path(trimmed_denovo), path(paf)
    output:
        path samplesheet_csv, emit: samplesheet_csv_de_novo_ch
        path trimmed_denovo, emit: trimmed_de_novo_fa_ch
        path paf, emit: paf_file_de_novo_ch
        path genbank_path, emit: genbank_path_ch
    script:
        samplesheet_csv = "entries.csv"
        // order is important, see Python script:
        """
        echo "$params.projectname,$entry,$barcode,$sample,$result,$genbank_path,$sample_fasta,$trimmed_denovo,$paf" >> $samplesheet_csv
        """    
}

process runReview_de_novo {
    publishDir 'results/dir3_review/n5_de_novo_results', mode: 'copy'

    input:
        file paf
        file genbank
        path trimmed_denovo
        path samplesheet_csv
    output:
        tuple path(pdf_file), path(samplesheet_csv)
    script:
        pdf_file = "de_novo_review.pdf"
        """
        review.py "$samplesheet_csv" "$params.plan_path_denovo" "$params.projectname" "$pdf_file"
        """
}


// Workflows:

workflow review_denovo {
    take: entries_de_novo_ch
    main:
        params.parts_path_denovo = file(params.all_parts)
        if( params.assembly_plan == "noplan" ) {
            params.plan_path_denovo = "noplan"
        }
        else {
            params.plan_path_denovo = file(params.assembly_plan)
        }
        convertGenbank_de_novo(entries_de_novo_ch)
        assembleDeNovo(convertGenbank_de_novo.out)
        trimAssembly(assembleDeNovo.out)
        alignParts_de_novo(trimAssembly.out, params.parts_path_denovo)
        writeCSV_de_novo(alignParts_de_novo.out)
        runReview_de_novo(writeCSV_de_novo.out.paf_file_de_novo_ch.collect(), writeCSV_de_novo.out.genbank_path_ch.collect(), writeCSV_de_novo.out.trimmed_de_novo_fa_ch.collect(), writeCSV_de_novo.out.samplesheet_csv_de_novo_ch.collectFile())
}
