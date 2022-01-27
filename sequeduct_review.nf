#!/usr/bin/env nextflow

params.results_csv = ''  // CSV of the sample ~ barcode relations. Columns: Barcode,Sample,Result,Review_consensus (may have other)
params.reference_dir = 'ref'  // dir of reference sequence Genbank files. Filenames (without extension) must match 'Sample' column entries

params.assembly_plan = ''  // Optional: the assembly plan CSV of the DNA constructs, one per line: Sample,Part_1,Part_2, etc. Must have a header line.
params.all_parts = ''  // FASTA file that contains all sequences to compare against


// params.fastq_filtered_dir = 'results/dir2_analysis/n2_fastq_filtered'  // The directory that contains the filtered FASTQ files
params.consensus_dir = 'results/dir2_analysis/n6_consensus'  // contains the consensus FASTA sequences

params.barcode_prefix = 'barcode'  // The prefix for the individual barcode directories as output by the sequencer

params.projectname = 'Noname'  // display in PDF

params.consensus_columname = 'Review_consensus'
params.consensus_true = 1  // marker for performing review


Channel
    .fromPath(params.results_csv)
    .splitCsv(header: true)
    // .unique { row -> row['Barcode_dir'] }
    .filter { item -> item['Review_consensus'] == params.consensus_true}
    
    .map { row -> 
        def barcode = row['Barcode']  // a barcode is present only once -- no pooling
        def sample = row['Sample']  // a sample may be present multiple times, in different barcodes
        def entry = "${barcode}_${sample}"  // unique key for each sample sheet entry
        def result = row["Result"]
        return [entry, barcode, sample, result]
        }
    .set { entries_ch }


process convertGenbank {

    input:
        tuple val(entry), val(barcode), val(sample), val(result) from entries_ch

    output:
        tuple val(entry), val(barcode), val(sample), val(result), path(sample_fasta) into entries_fasta_ch

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
        """
}
