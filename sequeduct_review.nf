#!/usr/bin/env nextflow

params.results_csv = ''  // CSV of the sample ~ barcode relations. Columns: Barcode,Sample,Result,Review_consensus (may have other)
params.reference_dir = 'ref'  // dir of reference sequence Genbank files. Filenames (without extension) must match 'Sample' column entries

params.assembly_plan = ''  // Optional: the assembly plan CSV of the DNA constructs, one per line: Sample,Part_1,Part_2, etc. Must have a header line.
params.all_parts = ''  // FASTA file that contains all sequences to compare against


params.fastq_filtered_dir = 'results/dir2_analysis/n2_fastq_filtered'  // The directory that contains the filtered FASTQ files
params.consensus_dir = 'results/dir2_analysis/n6_consensus'  // contains the consensus FASTA sequences

params.barcode_prefix = 'barcode'  // The prefix for the individual barcode directories as output by the sequencer

params.projectname = 'Noname'  // display in PDF

params.consensus_columname = 'Review_consensus'
params.consensus_true = '1'  // marker for performing review

params.denovo_columname = 'Review_de_novo'
params.denovo_true = '1'  // marker for performing review

assembly_prefix = 'egf'
canu_postfix = '.contigs.fasta'  // hardcoded into canu

///////////////////////////////////////////////////////////////////////////////
// Variant call consensus sequence review

Channel
    .fromPath(params.results_csv)
    .splitCsv(header: true)
    // .unique { row -> row['Barcode_dir'] }
    .filter { item -> item[params.consensus_columname] == params.consensus_true}
    
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
        tuple val(entry), val(barcode), val(sample), val(result), val(genbank_path), path(sample_fasta) into entries_fasta_ch

    script:
        genbank_path = PWD + '/' + params.reference_dir + '/' + sample + '.gb'
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
        """
}

process alignParts {
    publishDir 'results/dir3_review/n1_alignment', mode: 'copy', pattern: '*.paf'

    input:
        tuple val(entry), val(barcode), val(sample), val(result), val(genbank_path), path(sample_fasta) from entries_fasta_ch

    output:
        tuple val(entry), val(barcode), val(sample), val(result), val(genbank_path), path(sample_fasta), val(consensus_path), path(paf) into alignment_ch

    script:
        consensus_path = PWD + '/' + params.consensus_dir + '/' + entry + '_consensus.fa'
        parts_path = PWD + '/' + params.all_parts
        paf = entry + '.paf'
        """
        cat $sample_fasta $parts_path | \
        minimap2 -cx asm5 $consensus_path - > $paf
        """
}

process writeCSV {
    input:
        tuple val(entry), val(barcode), val(sample), val(result), val(genbank_path), path(sample_fasta), val(consensus_path), path(paf) from alignment_ch
    output:
        path samplesheet_csv into samplesheet_csv_ch
        path paf into paf_file_ch
    script:
        samplesheet_csv = "entries.csv"
        // order is important, see Python script:
        """
        echo "$params.projectname,$entry,$barcode,$sample,$result,$genbank_path,$sample_fasta,$consensus_path,$paf" >> $samplesheet_csv
        """    
}

process runReview {
    publishDir 'results/dir3_review/n2_results', mode: 'symlink'

    input:
        file paf from paf_file_ch.collect()
        path samplesheet_csv from samplesheet_csv_ch.collectFile()
    output:
        tuple path(pdf_file), path(samplesheet_csv), path(paf) into results_ch
    script:
        plan_path = PWD + '/' + params.assembly_plan
        pdf_file = "consensus_review.pdf"
        """
        #!/usr/bin/env python

        import os
        import pandas as pd
        import ediacara as edi

        entries = pd.read_csv("$samplesheet_csv", header=None)
        # see process writeCSV for columns:
        entries.columns = ['project', 'entry', 'barcode', 'sample', 'result', 'gb', 'fa', 'consensus', 'paf', ]

        consensus_list = []
        for index, row in entries.iterrows():
            assembly = edi.Assembly(assembly_path=row['consensus'],
                                    reference_path=row['gb'],
                                    alignment_path=row['paf'],
                                    assembly_plan="$plan_path")
            consensus_list += [assembly]

        assemblybatch = edi.AssemblyBatch(assemblies=consensus_list, name="$params.projectname")
        assemblybatch.perform_all_interpretations_in_group()

        edi.write_assembly_analysis_report("$pdf_file", assemblybatch)
        """
}


///////////////////////////////////////////////////////////////////////////////
// De novo assembly sequence review


Channel
    .fromPath(params.results_csv)
    .splitCsv(header: true)
    // .unique { row -> row['Barcode_dir'] }
    .filter { item -> item[params.denovo_columname] == params.denovo_true}
    
    .map { row -> 
        def barcode = row['Barcode']  // a barcode is present only once -- no pooling
        def sample = row['Sample']  // a sample may be present multiple times, in different barcodes
        def entry = "${barcode}_${sample}"  // unique key for each sample sheet entry
        def result = row["Result"]
        return [entry, barcode, sample, result]
        }
    .set { entries_de_novo_ch }

process convertGenbank_de_novo {
    input:
        tuple val(entry), val(barcode), val(sample), val(result) from entries_de_novo_ch

    output:
        tuple val(entry), val(barcode), val(sample), val(result), val(genbank_path), path(sample_fasta) into entries_fasta_de_novo_ch

    script:
        genbank_path = PWD + '/' + params.reference_dir + '/' + sample + '.gb'
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
        """
}

process assembleDeNovo {
    publishDir 'results/dir3_review/m1_de_novo_assembly', mode: 'symlink'
    
    input:
        tuple val(entry), val(barcode), val(sample), val(result), val(genbank_path), path(sample_fasta) from entries_fasta_de_novo_ch
    output:
        tuple val(entry), val(barcode), val(sample), val(result), val(genbank_path), path(sample_fasta), path(assembly_dir) into assembly_de_novo_ch 
    script:
        fastq_path = PWD + '/' + params.fastq_filtered_dir + '/' + barcode + '.fastq'
        assembly_dir = barcode + '_assembly'
        """
        canu -p $assembly_prefix -d $assembly_dir genomeSize=9k -nanopore $fastq_path
        """
}

process trimAssembly {
    publishDir 'results/dir3_review/o1_de_novo_assembly/trimmed', mode: 'symlink', pattern: '*_denovo.fasta'

    input:
        tuple val(entry), val(barcode), val(sample), val(result), val(genbank_path), path(sample_fasta), val(assembly_dir) from assembly_de_novo_ch 
    output:
        tuple val(entry), val(barcode), val(sample), val(result), val(genbank_path), path(sample_fasta), val(assembly_dir), path(trimmed_denovo) into trimmed_de_novo_ch 
    
    script:
        trimmed_denovo = barcode + '_denovo.fasta'
        """
        #!/usr/bin/env python
        from Bio import SeqIO

        canu_fasta = "$assembly_dir" + '/' + "$assembly_prefix" + "$canu_postfix"
        try:
            contig = SeqIO.read(canu_fasta, format="fasta")
        except:
            print("The FASTA file contains more than 1 contigs. First contig used.")
            contig = next(SeqIO.parse(canu_fasta, format="fasta"))

        entries = contig.description.split(" ")
        desc_dict = {"name": entries[0]}  # first is the name
        for entry in entries[1:]:  # addressed the first one above
            elements = entry.split("=")
            desc_dict[elements[0]] = elements[1]

        if desc_dict["suggestCircular"] == "yes":  # as output by canu
            start, end = desc_dict["trim"].split("-")  # must contain 2 values
            start = int(start)
            end = int(end)
            SeqIO.write(contig[start:end], "$trimmed_denovo", format="fasta")
        else:  # keep intact
            SeqIO.write(contig, "$trimmed_denovo", format="fasta")
        
        print("Trimmed:", "$barcode")
        """
}
