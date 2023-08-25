#!/usr/bin/env python
# Copyright 2021 Edinburgh Genome Foundry, University of Edinburgh
#
# This file is part of Sequeduct.
#
# Sequeduct is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
#
# Sequeduct is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License along with Sequeduct. If not, see <https:www.gnu.org/licenses/>.

import resource
import sys

# Set 5 times the default. This is to avoid a problem during PDF build and may be revised later.
resource.setrlimit(resource.RLIMIT_STACK, (41943040, -1))
sys.setrecursionlimit(5000)

print("RECURSION LIMIT:")
print(sys.getrecursionlimit())
print()
print()

samplesheet_csv = sys.argv[1]  # skip first filename
params_projectname = sys.argv[2]
pdf_file = sys.argv[3]
results_csv_file = sys.argv[4]
low_depth_value = sys.argv[5]

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
sequencinggroup = edi.SequencingGroup(
    comparatorgroups, name=params_projectname, low_depth_cutoff=low_depth_value
)
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
