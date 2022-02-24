#!/usr/bin/env nextflow
nextflow.enable.dsl=2

params.sample_sheet = ''  // CSV of the sample ~ barcode relations. Columns: Sample,Barcode_dir (may have other)
params.fastq_dir = 'fastq'  // The directory that contains the barcode directories of FASTQ files
params.barcode_prefix = 'barcode'  // The prefix for the individual barcode directories as output by the sequencer

params.reference_dir = ''  // dir of reference sequence Genbank files. Filenames (without extension) must match 'Sample' column entries
params.projectname = 'Noname'  // display in PDF


//
params.max_len_fraction = 1.5  // For calculating max length of filtered FASTQ reads, to avoid plasmid dimers
params.quality_cutoff = 9  // For NanoFilt
params.min_length = 500  // nucleotides
//

params.assembly_plan = ''  // Optional: the assembly plan CSV of the DNA constructs, one per line: Sample,Part_1,Part_2, etc. Must have a header line.
params.all_parts = ''  // FASTA file that contains all sequences to compare against

params.fastq_filtered_dir = 'results/dir2_analysis/n2_fastq_filtered'  // The directory that contains the filtered FASTQ files
params.consensus_dir = 'results/dir2_analysis/n6_consensus'  // contains the consensus FASTA sequences

params.consensus_columname = 'Review_consensus'
params.consensus_true = '1'  // marker for performing review

params.denovo_columname = 'Review_de_novo'
params.denovo_true = '1'  // marker for performing review

params.assembly_prefix = 'egf'
params.canu_postfix = '.contigs.fasta'  // hardcoded into canu

//

include { preview_workflow } from "$projectDir/nextflow/sequeduct_preview.nf"
include { analysis_workflow } from "$projectDir/nextflow/sequeduct_analysis.nf"
include { review_consensus } from "$projectDir/nextflow/sequeduct_review.nf"
include { review_denovo } from "$projectDir/nextflow/sequeduct_review.nf"



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
