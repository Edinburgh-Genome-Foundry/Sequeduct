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

genbank_path = sys.argv[1]  # skip first filename
sample = sys.argv[2]
sample_fasta = sys.argv[3]
max_len_fraction = sys.argv[4]
flag = sys.argv[5]

from Bio import SeqIO

# Genbank in
record = SeqIO.read(genbank_path, "genbank")
record.id = sample
# FASTA out
with open(sample_fasta, "w") as output_handle:
    SeqIO.write(record, output_handle, "fasta")

if flag == "none":
    pass  # no need for stdout
elif flag == "length":
    cutoff_length = int(float(len(record)) * float(max_len_fraction))
    print(cutoff_length, end="")
elif flag == "canu":
    print(str(round(len(record) / 1000)), end="")  # to get k value for canu
