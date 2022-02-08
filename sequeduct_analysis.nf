#!/usr/bin/env nextflow

params.sample_sheet = ''  // CSV of the sample ~ barcode relations. Columns: Sample,Barcode_dir (may have other)
params.reference_dir = ''  // dir of reference sequence Genbank files. Filenames (without extension) must match 'Sample' column entries

params.fastq_dir = 'fastq'  // The directory that contains the barcode directories of FASTQ files
params.barcode_prefix = 'barcode'  // The prefix for the individual barcode directories as output by the sequencer
params.max_len_fraction = 1.5  // For calculating max length of filtered FASTQ reads, to avoid plasmid dimers

params.projectname = 'Noname'  // display in PDF

quality_cutoff = 9
min_length = 500  // nucleotides

ref_dir = file(params.reference_dir)


Channel
    .fromPath(params.sample_sheet)
    .splitCsv(header: true)
    // .unique { row -> row['Barcode_dir'] }
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


process convertGenbank {
    publishDir 'results/dir2_analysis/n1_ref_fasta', mode: 'copy', pattern: '*.fa'  // need only the fasta

    input:
        tuple val(entry), val(barcode), file(barcode_path), val(fastq_files), val(sample), file(genbank_path) from entries_ch

    output:
        tuple val(entry), val(barcode), file(barcode_path), val(fastq_files), val(sample), path(sample_fasta), stdout into entries_fasta_ch  // stdout for seq length

    script:
        sample_fasta = sample + '.fa'

        """
        #!/usr/bin/env python

        import os
        from Bio import SeqIO

        # Genbank in
        record = SeqIO.read("$genbank_path", "genbank")
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
        tuple val(entry), val(barcode), path(barcode_path), val(fastq_files), val(sample), path(sample_fasta), val(seq_length) from entries_fasta_ch

    output:
        tuple val(entry), val(barcode), path(fastq_file) into fastq_filtered_ch
        tuple val(entry), val(barcode), val(sample), path(sample_fasta), val(seq_length), path(fastq_file) into entries_fasta_fastq_ch

    script:
        fastq_file = barcode + '.fastq'  // need for output

        fastqFileString = fastq_files.join(' ')  // need as one string for cat

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
        tuple val(entry), val(barcode), val(sample), path(sample_fasta), val(seq_length), path(fastq_file), path(paf_file), path(bam_file), path(bai_file), path(counts_tsv) into entries_aligned_ch

    script:
        sam_file = entry + '.sam'
        paf_file = entry + '.paf'
        sorted_sam_file = entry + '_sorted.sam'
        counts_tsv = entry + '_sorted.sam_counts.tsv'
        bam_file = entry + '_sorted.bam'
        bai_file = entry + '_sorted.bam.bai'  // index file
        """
        minimap2 -a $sample_fasta $fastq_file > $sam_file
        paftools.js sam2paf $sam_file > $paf_file
        samtools sort -O sam -T sample.sort -o $sorted_sam_file $sam_file
        samtools depth -aa $sorted_sam_file > $counts_tsv
        samtools view -S -b $sorted_sam_file > $bam_file
        samtools index $bam_file $bai_file
        """
}


process callVariants {
    publishDir 'results/dir2_analysis/n5_variant_calls', mode: 'copy', pattern: '*.vcf'

    input:
        tuple val(entry), val(barcode), val(sample), path(sample_fasta), val(seq_length), path(fastq_file), path(paf_file), path(bam_file), path(bai_file), path(counts_tsv) from entries_aligned_ch

    output:
        tuple val(entry), val(barcode), val(sample), path(sample_fasta), val(seq_length), path(fastq_file), path(paf_file), path(bam_file), path(bai_file), path(counts_tsv), path(vcf_file) into entries_vcf_ch

    script:
        vcf_file = entry + '.vcf'
        """
        freebayes --ploidy 1 --min-alternate-fraction 0.1 --min-alternate-count 2 -f $sample_fasta $bam_file > $vcf_file
        """
}

process callConsensus {
    publishDir 'results/dir2_analysis/n6_consensus', mode: 'copy', pattern: '*_consensus.fa'

    input:
        tuple val(entry), val(barcode), val(sample), path(sample_fasta), val(seq_length), path(fastq_file), path(paf_file), path(bam_file), path(bai_file), path(counts_tsv), path(vcf_file) from entries_vcf_ch

    output:
        tuple val(entry), val(barcode), val(sample), path(sample_fasta), val(seq_length), path(fastq_file), path(paf_file), path(bam_file), path(bai_file), path(counts_tsv), path(vcf_file), path(filtered_vcf_file), path(consensus_fa_file) into entries_out_ch

    script:
        vcf_gz_file = entry + '.vcf.gz'
        filtered_vcf_file = entry + '_filtered.vcf'
        filtered_vcf_gz_file = entry + '_filtered.vcf.gz'
        consensus_fa_file = entry + '_consensus.fa'
        """
        bgzip --keep --index $vcf_file
        bcftools index $vcf_gz_file
        bcftools filter --output-type v -i'%QUAL>10' $vcf_gz_file > $filtered_vcf_file
        bgzip --keep --index $filtered_vcf_file
        bcftools index $filtered_vcf_gz_file
        bcftools consensus --fasta-ref $sample_fasta --output $consensus_fa_file $filtered_vcf_gz_file
        """
}

process writeCSV {
    input:
        tuple val(entry), val(barcode), val(sample), path(sample_fasta), val(seq_length), path(fastq_file), path(paf_file), path(bam_file), path(bai_file), path(counts_tsv), path(vcf_file), path(filtered_vcf_file), path(consensus_fa_file) from entries_out_ch
    output:
        path samplesheet_csv into samplesheet_csv_ch
        path paf_file into paf_file_ch
        path counts_tsv into counts_tsv_ch
        path filtered_vcf_file into filtered_vcf_file_ch
        path consensus_fa_file into consensus_fa_file_ch
    script:
        samplesheet_csv = "entries.csv"
        // order is important, see Python script:
        """
        echo "$params.projectname,$entry,$barcode,$sample,$sample_fasta,$filtered_vcf_file,$paf_file,$counts_tsv,$consensus_fa_file" >> $samplesheet_csv
        """    
}


Channel
    .fromPath(params.sample_sheet)
    .splitCsv(header: true)
    .unique { row -> row['Sample'] }
    .map { row -> file(params.reference_dir + '/' + row['Sample'] + '.gb') }
    .set { genbank_ch }

process runEdiacara {
    publishDir 'results/dir2_analysis/n7_results', mode: 'copy'

    input:
        file paf from paf_file_ch.collect()
        file counts_tsv from counts_tsv_ch.collect()
        file filtered_vcf_file from filtered_vcf_file_ch.collect()
        file consensus_fa_file from consensus_fa_file_ch.collect()
        path genbank from genbank_ch.collect()
        path samplesheet_csv from samplesheet_csv_ch.collectFile()
    output:
        tuple path(pdf_file), path(results_csv_file), path(samplesheet_csv) into results_ch
    script:
        pdf_file = "Ediacara_report.pdf"
        results_csv_file = "results.csv"
        """
        #!/usr/bin/env python

        import os
        import pandas as pd
        from Bio import SeqIO
        import ediacara as edi

        entries = pd.read_csv('$samplesheet_csv', header=None)
        # match writeCSV process:
        entries.columns = ['projectname', 'entry', 'barcode', 'sample', 'fasta', 'vcf', 'paf', 'tsv', 'consensus_fasta']

        comparatorgroups = []
        for index, row in entries.iterrows():
            print('Processing', row['entry'], end='')
            entry = row['entry']
            sample = row['sample']
            vcf = row['vcf']

            reference_gb = row['sample'] + '.gb'
            record = SeqIO.read(reference_gb, 'genbank')
            record.id = sample
            references = {record.id: record}

            tsv_file = row['tsv']
            paf_path = row['paf']

            tsv = edi.ComparatorGroup.load_tsv(tsv_file)
            paf = edi.ComparatorGroup.load_paf(paf_path)

            assembly_paths = {
                sample: row['consensus_fasta']
            }
            vcf_paths = {
                sample: vcf
            }

            comparator_group = edi.ComparatorGroup(references=references,
                                                alignments={'paf': paf, 'tsv': tsv},
                                                barcode=row['barcode'],
                                                assembly_paths=assembly_paths,
                                                vcf_paths=vcf_paths)

            list_of_constructs = [sample]
            for element in list_of_constructs:
                comparator_group.add_comparator(element)

            comparatorgroups += [comparator_group]
            print('    ... done')


        # Create PDF report
        sequencinggroup = edi.SequencingGroup(comparatorgroups, name='$params.projectname')
        sequencinggroup.perform_all_comparisons_in_sequencinggroup()
        edi.write_sequencinggroup_report(target='$pdf_file', sequencinggroup=sequencinggroup)

        print('PDF created')
        ###############################################################################
        barcodes = []
        samples = []
        results = []

        for comparatorgroup in sequencinggroup.comparatorgroups:
            for index, row in comparatorgroup.summary_table.iterrows():
                barcodes.append(comparatorgroup.barcode)
                samples += [row["Name"]]
                results += [row["Result"]]

        d = {
            'Barcode': pd.Series(barcodes),
            'Sample': pd.Series(samples),
            'Result': pd.Series(results),
        }

        results_table = pd.DataFrame(d)

        results_table['Review_consensus'] = 0
        results_table['Review_de_novo'] = 0

        results_table.to_csv('$results_csv_file', index=False)
        print('Done')
        """
}
