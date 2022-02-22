#!/usr/bin/env nextflow
nextflow.enable.dsl=2

params.results_csv = ''  // CSV of the sample ~ barcode relations. Columns: Barcode,Sample,Result,Review_consensus,Review_de_novo (may have other)
params.reference_dir = ''  // dir of reference sequence Genbank files. Filenames (without extension) must match 'Sample' column entries

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


parts_path = file(params.all_parts)
plan_path = file(params.assembly_plan)

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
        #!/usr/bin/env python

        import os
        import pandas as pd
        import ediacara as edi

        entries = pd.read_csv("$samplesheet_csv", header=None)
        # see process writeCSV for columns:
        entries.columns = ['project', 'entry', 'barcode', 'sample', 'result', 'gb', 'fa', 'consensus', 'paf', ]
        entries.sort_values(by=['barcode', 'sample'], inplace=True)  # have them in order in the pdf

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


workflow review_consensus {
    take: entries_ch
    main:
        convertGenbank(entries_ch)
        alignParts(convertGenbank.out, parts_path)
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
        #!/usr/bin/env python

        import os
        from Bio import SeqIO

        # Genbank in
        record = SeqIO.read("$genbank_path", "genbank")
        record.id = "$sample"
        # FASTA out
        with open("$sample_fasta", "w") as output_handle:
            SeqIO.write(record, output_handle, "fasta")

        print(str(round(len(record) / 1000)), end='')  # to get k value for canu
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
        canu -p $assembly_prefix -d $assembly_dir $genomsize_param -nanopore $fastq_path
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
        #!/usr/bin/env python

        import os
        import pandas as pd
        import ediacara as edi

        entries = pd.read_csv("$samplesheet_csv", header=None)
        # see process writeCSV for columns:
        entries.columns = ['project', 'entry', 'barcode', 'sample', 'result', 'gb', 'fa', 'de_novo', 'paf', ]
        entries.sort_values(by=['barcode', 'sample'], inplace=True)  # have them in order in the pdf

        consensus_list = []
        for index, row in entries.iterrows():
            assembly = edi.Assembly(assembly_path=row['de_novo'],
                                    reference_path=row['gb'],
                                    alignment_path=row['paf'],
                                    assembly_plan="$plan_path")
            consensus_list += [assembly]

        assemblybatch = edi.AssemblyBatch(assemblies=consensus_list, name="$params.projectname")
        assemblybatch.perform_all_interpretations_in_group()

        edi.write_assembly_analysis_report("$pdf_file", assemblybatch)
        """
}


workflow review_denovo {
    take: entries_de_novo_ch
    main:
        convertGenbank_de_novo(entries_de_novo_ch)
        assembleDeNovo(convertGenbank_de_novo.out)
        trimAssembly(assembleDeNovo.out)
        alignParts_de_novo(trimAssembly.out, parts_path)
        writeCSV_de_novo(alignParts_de_novo.out)
        runReview_de_novo(writeCSV_de_novo.out.paf_file_de_novo_ch.collect(), writeCSV_de_novo.out.genbank_path_ch.collect(), writeCSV_de_novo.out.trimmed_de_novo_fa_ch.collect(), writeCSV_de_novo.out.samplesheet_csv_de_novo_ch.collectFile())
}


workflow {

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
        // .unique { row -> row['Barcode_dir'] }
        .filter { item -> item[params.denovo_columname] == params.denovo_true}
        
        .map { row -> 
            def barcode = row['Barcode']  // a barcode is present only once -- no pooling
            def sample = row['Sample']  // a sample may be present multiple times, in different barcodes
            def entry = "${barcode}_${sample}"  // unique key for each sample sheet entry
            def result = row["Result"]
            def genbank_path = file("${params.reference_dir}/${sample}.gb")
            def fastq_path = file("${params.fastq_filtered_dir}/${barcode}.fastq")
            // def assembly_dir = file("${params.reference_dir}/${sample}.gb")
            return [entry, barcode, sample, result, genbank_path, fastq_path]
            }
        .set { entries_de_novo_ch }


    review_consensus(entries_ch)

    review_denovo(entries_de_novo_ch)
}