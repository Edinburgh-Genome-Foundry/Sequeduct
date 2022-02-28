#!/usr/bin/env nextflow
nextflow.enable.dsl=2

process convertGenbank {
    publishDir 'results/dir2_analysis/n1_ref_fasta', mode: 'copy', pattern: '*.fa'  // need only the fasta

    input:
        tuple val(entry), val(barcode), file(barcode_path), val(fastq_files), val(sample), file(genbank_path)

    output:
        tuple val(entry), val(barcode), file(barcode_path), val(fastq_files), val(sample), path(sample_fasta), stdout // stdout for seq length

    script:
        sample_fasta = sample + '.fa'

        """
        convert_genbank.py "$genbank_path" "$sample" "$sample_fasta" "length"
        """
}

process runNanoFilt {
    publishDir 'results/dir2_analysis/n2_fastq_filtered', mode: 'copy', pattern: '*.fastq'  // need only the fastq

    input:
        tuple val(entry), val(barcode), path(barcode_path), val(fastq_files), val(sample), path(sample_fasta), val(seq_length)

    output:
        tuple val(entry), val(barcode), path(fastq_file), emit: fastq_filtered_ch
        tuple val(entry), val(barcode), val(sample), path(sample_fasta), val(seq_length), path(fastq_file), emit: entries_fasta_fastq_ch

    script:
        fastq_file = barcode + '.fastq'  // need for output

        fastqFileString = fastq_files.join(' ')  // need as one string for cat

        max_length = seq_length * params.max_len_fraction
        """
        cat $fastqFileString | NanoFilt -l $params.min_length --maxlength $max_length -q $params.quality_cutoff > $fastq_file
        """
}

process runNanoPlot {
    publishDir 'results/dir2_analysis/n3_nanoplots', mode: 'copy'

    input:
        tuple val(entry), val(barcode), path(fastq_file)

    output:
        path barcode

    script:
        """
        NanoPlot --raw --fastq $fastq_file -o $barcode
        """
}

process alignEntries {
    publishDir 'results/dir2_analysis/n4_alignment', mode: 'copy', pattern: '*.paf'
    publishDir 'results/dir2_analysis/n4_alignment', mode: 'copy', pattern: '*.bam'
    publishDir 'results/dir2_analysis/n4_alignment', mode: 'copy', pattern: '*.bai'
    publishDir 'results/dir2_analysis/n4_alignment', mode: 'copy', pattern: '*.tsv'

    input:
        tuple val(entry), val(barcode), val(sample), path(sample_fasta), val(seq_length), path(fastq_file)

    output:
        tuple val(entry), val(barcode), val(sample), path(sample_fasta), val(seq_length), path(fastq_file), path(paf_file), path(bam_file), path(bai_file), path(counts_tsv)

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
        tuple val(entry), val(barcode), val(sample), path(sample_fasta), val(seq_length), path(fastq_file), path(paf_file), path(bam_file), path(bai_file), path(counts_tsv)

    output:
        tuple val(entry), val(barcode), val(sample), path(sample_fasta), val(seq_length), path(fastq_file), path(paf_file), path(bam_file), path(bai_file), path(counts_tsv), path(vcf_file)

    script:
        vcf_file = entry + '.vcf'
        """
        freebayes --ploidy 1 --min-alternate-fraction 0.1 --min-alternate-count 2 -f $sample_fasta $bam_file > $vcf_file
        """
}

process callConsensus {
    publishDir 'results/dir2_analysis/n6_consensus', mode: 'copy', pattern: '*_consensus.fa'

    input:
        tuple val(entry), val(barcode), val(sample), path(sample_fasta), val(seq_length), path(fastq_file), path(paf_file), path(bam_file), path(bai_file), path(counts_tsv), path(vcf_file)

    output:
        tuple val(entry), val(barcode), val(sample), path(sample_fasta), val(seq_length), path(fastq_file), path(paf_file), path(bam_file), path(bai_file), path(counts_tsv), path(vcf_file), path(filtered_vcf_file), path(consensus_fa_file)

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
        tuple val(entry), val(barcode), val(sample), path(sample_fasta), val(seq_length), path(fastq_file), path(paf_file), path(bam_file), path(bai_file), path(counts_tsv), path(vcf_file), path(filtered_vcf_file), path(consensus_fa_file)
    output:
        path samplesheet_csv, emit: samplesheet_csv_ch
        path paf_file, emit: paf_file_ch
        path counts_tsv, emit: counts_tsv_ch
        path filtered_vcf_file, emit: filtered_vcf_file_ch
        path consensus_fa_file, emit: consensus_fa_file_ch
    script:
        samplesheet_csv = "entries.csv"
        // order is important, see Python script:
        """
        echo "$params.projectname,$entry,$barcode,$sample,$sample_fasta,$filtered_vcf_file,$paf_file,$counts_tsv,$consensus_fa_file" >> $samplesheet_csv
        """    
}

process runEdiacara {
    publishDir 'results/dir2_analysis/n7_results', mode: 'copy'

    input:
        file paf
        file counts_tsv
        file filtered_vcf_file
        file consensus_fa_file
        path genbank
        path samplesheet_csv
    output:
        tuple path(pdf_file), path(results_csv_file), path(samplesheet_csv)
    script:
        pdf_file = "Ediacara_report.pdf"
        results_csv_file = "results.csv"
        """
        analysis.py $samplesheet_csv "$params.projectname" $pdf_file $results_csv_file
        """
}


workflow analysis_workflow {
    take:
        entries_ch
        genbank_ch
    main:
        convertGenbank(entries_ch)
        runNanoFilt(convertGenbank.out)
        runNanoPlot(runNanoFilt.out.fastq_filtered_ch)
        alignEntries(runNanoFilt.out.entries_fasta_fastq_ch)
        callVariants(alignEntries.out)
        callConsensus(callVariants.out)
        writeCSV(callConsensus.out)
        runEdiacara(writeCSV.out.paf_file_ch.collect(), writeCSV.out.counts_tsv_ch.collect(), writeCSV.out.filtered_vcf_file_ch.collect(), writeCSV.out.consensus_fa_file_ch.collect(), genbank_ch.collect(), writeCSV.out.samplesheet_csv_ch.collectFile())
}
