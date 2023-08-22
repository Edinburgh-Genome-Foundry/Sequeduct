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
