#!/usr/bin/env python
import sys

samplesheet_csv = sys.argv[1]  # skip first filename
params_projectname = sys.argv[2]
pdf_file = sys.argv[3]
results_csv_file = sys.argv[4]

import pandas as pd
from Bio import SeqIO
import ediacara as edi

entries = pd.read_csv(samplesheet_csv, header=None)
# match writeCSV process:
entries.columns = [
    "projectname",
    "entry",
    "barcode",
    "sample",
    "fasta",
    "vcf",
    "paf",
    "tsv",
    "consensus_fasta",
]
entries.sort_values(
    by=["barcode", "sample"], inplace=True
)  # have them in order in the pdf

comparatorgroups = []
for index, row in entries.iterrows():
    print("Processing", row["entry"], end="")
    entry = row["entry"]
    sample = row["sample"]
    vcf = row["vcf"]

    reference_gb = row["sample"] + ".gb"
    record = SeqIO.read(reference_gb, "genbank")
    record.id = sample
    references = {record.id: record}

    tsv_file = row["tsv"]
    paf_path = row["paf"]

    tsv = edi.ComparatorGroup.load_tsv(tsv_file)
    paf = edi.ComparatorGroup.load_paf(paf_path)

    assembly_paths = {sample: row["consensus_fasta"]}
    vcf_paths = {sample: vcf}

    comparator_group = edi.ComparatorGroup(
        references=references,
        alignments={"paf": paf, "tsv": tsv},
        barcode=row["barcode"],
        assembly_paths=assembly_paths,
        vcf_paths=vcf_paths,
    )

    list_of_constructs = [sample]
    for element in list_of_constructs:
        comparator_group.add_comparator(element)

    comparatorgroups += [comparator_group]
    print("    ... done")


# Create PDF report
sequencinggroup = edi.SequencingGroup(comparatorgroups, name=params_projectname)
sequencinggroup.perform_all_comparisons_in_sequencinggroup()
edi.write_sequencinggroup_report(target=pdf_file, sequencinggroup=sequencinggroup)

print("PDF created")
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
    "Barcode": pd.Series(barcodes),
    "Sample": pd.Series(samples),
    "Result": pd.Series(results),
}

results_table = pd.DataFrame(d)

results_table["Review_consensus"] = 0
results_table["Review_de_novo"] = 0

results_table.to_csv(results_csv_file, index=False)
print("Done")
