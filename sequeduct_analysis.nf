#!/usr/bin/env nextflow

params.sample_sheet = ''  // CSV of the sample ~ barcode relations. Columns: Sample,Barcode_dir (may have other)
params.reference_dir = 'ref'   // dir of reference sequence Genbank files. Filenames must match 'Sample' column entries

params.fastq_dir = 'fastq'  // The directory that contains the barcode directories of FASTQ files
params.barcode_prefix = 'barcode'  // The prefix for the individual barcode directories as output by the sequencer
params.max_len_fraction = 1.5  // For calculating max length of filtered FASTQ reads, to avoid plasmid dimers

quality_cutoff = 9
min_length = 500  // nucleotides


Channel
    .fromPath(params.sample_sheet)
    .splitCsv(header: true)
    // .unique { row -> row['Barcode_dir'] }
    .map { row -> 
        def barcode_dir = row['Barcode_dir']  // a barcode is present only once -- no pooling
        def sample = row['Sample']  // a sample may be present multiple times, in different barcodes
        def entry = "${barcode_dir}_${sample}"  // unique key for each sample sheet entry
        return [entry, barcode_dir, sample]
        }
    .set { entries_ch  }


process convertGenbank {
    publishDir 'results/dir2_analysis/n1_ref_fasta', mode: 'copy', pattern: '*.fa'  // need only the fasta

    input:
        tuple val(entry), val(barcode), val(sample) from entries_ch

    output:
        tuple val(entry), val(barcode), val(sample), path(sample_fasta), stdout into entries_fasta_ch  // stdout for seq length

    script:
        genbank_path = params.reference_dir + '/' + sample + '.gb'
        sample_fasta = sample + '.fa'

        """
        #!/usr/bin/env python

        import os
        from Bio import SeqIO

        # Genbank in
        record = SeqIO.read(os.path.join("$PWD", "$genbank_path"), "genbank")
        record.id = "$sample"
        # FASTA out
        with open("$sample_fasta", "w") as output_handle:
            SeqIO.write(record, output_handle, "fasta")

        print(len(record), end="")
        """
}


process runNanoFilt {
    publishDir 'results/dir2_analysis/n2_fastq_filtered', mode: 'copy', pattern: '*.fastq'  // need only the fastq

    input:
        tuple val(entry), val(barcode), val(sample), path(sample_fasta), val(seq_length) from entries_fasta_ch

    output:
        tuple val(entry), val(barcode), path(fastq_file) into fastq_filtered_ch
        tuple val(entry), val(barcode), val(sample), path(sample_fasta), val(seq_length), path(fastq_file) into entries_fasta_fastq_ch

    script:
        barcode_path = params.fastq_dir + '/' + barcode
        fastq_file = barcode + '.fastq'  // need for output
        fastqDir = file(barcode_path)

        fastqFiles = fastqDir.listFiles()  // multiple FASTQ in each barcode
        fastqFilePaths = []
        fastqFileString = fastqFiles.join(' ')  // need as one string for cat

        max_length = seq_length * params.max_len_fraction
        """
        cat $fastqFileString | NanoFilt -l $min_length --maxlength $max_length -q $quality_cutoff > $fastq_file
        """
}

process runNanoPlot {
    publishDir 'results/dir2_analysis/n3_nanoplots', mode: 'copy'

    input:
        tuple val(entry), val(barcode), path(fastq_file) from fastq_filtered_ch

    output:
        path barcode into nanoplots
        stdout result

    script:
        """
        NanoPlot --raw --fastq $fastq_file -o $barcode
        """
}

process AlignEntries {
    publishDir 'results/dir2_analysis/n4_alignment', mode: 'copy'

    input: 
        tuple val(entry), val(barcode), val(sample), path(sample_fasta), val(seq_length), path(fastq_file) from entries_fasta_fastq_ch

    output:
        path sam_file
        path paf_file
        path sorted_sam_file

    script:
        sam_file = entry + '.sam'
        paf_file = entry + '.paf'
        sorted_sam_file = entry + '_sorted.sam'
        """
        minimap2 -a $sample_fasta $fastq_file > $sam_file
        paftools.js sam2paf $sam_file > $paf_file
        samtools sort -O sam -T sample.sort -o $sorted_sam_file $sam_file
        """
}

result.view()
