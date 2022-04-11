#!/usr/bin/env python
import sys

samplesheet_csv = sys.argv[1]  # skip first filename
params_plan_path = sys.argv[2]
params_projectname = sys.argv[3]
pdf_file = sys.argv[4]

import pandas as pd
import ediacara as edi

entries = pd.read_csv(samplesheet_csv, header=None)
# see process writeCSV for columns:
entries.columns = [
    "project",
    "entry",
    "barcode",
    "sample",
    "result",
    "gb",
    "fa",
    "consensus",
    "paf",
]
entries.sort_values(
    by=["barcode", "sample"], inplace=True
)  # have them in order in the pdf

if params_plan_path == "noplan":
    assembly_plan = None
else:
    assembly_plan = params_plan_path
consensus_list = []
for index, row in entries.iterrows():
    assembly = edi.Assembly(
        assembly_path=row["consensus"],
        reference_path=row["gb"],
        alignment_path=row["paf"],
        assembly_plan=assembly_plan,
    )
    consensus_list += [assembly]

assemblybatch = edi.AssemblyBatch(assemblies=consensus_list, name=params_projectname)
assemblybatch.perform_all_interpretations_in_group()

edi.write_assembly_analysis_report(pdf_file, assemblybatch)
